function runChunkedPopavgPerBehavior_plotOnly_getMouseDataNames_FIXED( ...
    mouseIDs, baseSessionNames, probeRegions, saveFile)
% Corrected plot-only companion for the combined CHUNKED trial-averaged
% population-average behavior-specific analysis.
%
% IMPORTANT FIXES:
%   1) Animals are explicitly mapped from the requested mouseIDs to the
%      saved allSessions.mouseIDs before ANY values are extracted.
%      This prevents a saved-order/requested-order mismatch from assigning
%      one animal's real lag or permutation CI to another animal.
%
%   2) The 95% permutation interval is calculated ONCE and stored in
%      lagCILowMat / lagCIHighMat. The dot-plot interval lines and the
%      printed table use these exact same stored numbers.
%
%   3) The CI plot draws explicit endpoint markers at the exact lower and
%      upper percentile values, so the plotted bounds can be visually
%      checked directly against PermCI_Low_s and PermCI_High_s in the table.
%
%   4) Accepted permutations are handled safely as:
%         accepted == 1 AND finite peak lag
%      rather than blindly converting arbitrary numeric values to logical.
%
% Expected combined structure:
%   allSessions.behaviorResults{bIdx}
%
% Each behavior result contains:
%   B.real_xc                 [nAnimals x nLags]
%   B.real_peakLagSec         [nAnimals x 1]
%   B.real_peakCorr           [nAnimals x 1]
%   B.real_nTrialsUsed        [nAnimals x 1]
%   B.perm_peakLagSec         [nAnimals x 300]
%   B.perm_accepted           [nAnimals x 300]
%   B.shiftCorrUpper95        [nAnimals x 1]
%
% Run:
%   runChunkedPopavgPerBehavior_plotOnly_getMouseDataNames_FIXED

%% ---------------- defaults ----------------

if nargin < 1 || isempty(mouseIDs)
    mouseIDs = { ...
        'D026', ...
        'D020', ...
        'D024', ...
        'D043', ...
        'D050', ...
        'D054'};
end

if nargin < 2 || isempty(baseSessionNames)
    baseSessionNames = { ...
        'D026-032923-ArenaRecording', ...
        'D020-062922-ArenaRecording', ...
        'D024-111022-ArenaRecording', ...
        'D043-020525-ArenaRecording', ...
        'D050-121825-ArenaRecording', ...
        'D054-012126-ArenaRecording'};
end

if nargin < 3 || isempty(probeRegions)
    probeRegions = { ...
        '', ...
        '', ...
        '', ...
        '', ...
        'CFA', ...
        'CFA'};
end

if nargin < 4 || isempty(saveFile)
    saveFile = fullfile( ...
        'C:\Users\mirilab\Documents\GlobusTransfer', ...
        'concatCrossCorrPerCanonicalBehavior_classifier_chunked_trialavg_popavg_ALL_SESSIONS.mat');
end

nSess = numel(mouseIDs);

if numel(baseSessionNames) ~= nSess || numel(probeRegions) ~= nSess
    error('mouseIDs, baseSessionNames, and probeRegions must have the same length.');
end

%% ---------------- use David's getMouseDataNames ----------------

resolvedDataNames = cell(1,nSess);

for s = 1:nSess
    resolvedDataNames{s} = getMouseDataNames( ...
        mouseIDs{s}, ...
        baseSessionNames{s}, ...
        probeRegions{s});
end

%% ---------------- load combined analysis ----------------

if ~isfile(saveFile)
    error('Saved combined chunked results file not found:\n%s',saveFile);
end

S = load(saveFile,'allSessions');

if ~isfield(S,'allSessions')
    error('File does not contain variable "allSessions":\n%s',saveFile);
end

allSessions = S.allSessions;

if ~isfield(allSessions,'behaviorResults') || isempty(allSessions.behaviorResults)
    error('allSessions.behaviorResults is missing or empty.');
end

behaviorResults = allSessions.behaviorResults;

if isfield(allSessions,'behaviors') && ~isempty(allSessions.behaviors)
    behaviors = allSessions.behaviors(:)';
else
    behaviors = nan(1,numel(behaviorResults));
    for bIdx = 1:numel(behaviorResults)
        behaviors(bIdx) = behaviorResults{bIdx}.beh;
    end
end

nBeh = numel(behaviors);

%% ========================================================================
%% FIX 1: explicitly map requested animal order -> saved animal row
%% ========================================================================

animalLabels = string(mouseIDs(:))';

if isfield(allSessions,'mouseIDs') && ~isempty(allSessions.mouseIDs)

    savedMouseIDs = string(allSessions.mouseIDs(:));
    savedRowForRequestedAnimal = nan(nSess,1);

    for s = 1:nSess

        hit = find(savedMouseIDs == string(mouseIDs{s}),1,'first');

        if isempty(hit)
            error(['Requested mouse %s was not found in allSessions.mouseIDs.\n' ...
                   'Saved mice are: %s'], ...
                   mouseIDs{s},strjoin(cellstr(savedMouseIDs),', '));
        end

        savedRowForRequestedAnimal(s) = hit;
    end

else
    warning(['allSessions.mouseIDs is unavailable. Assuming saved row order ' ...
             'matches supplied mouseIDs exactly.']);
    savedRowForRequestedAnimal = (1:nSess)';
end

fprintf('\nSaved-row mapping used for plotting:\n');
for s = 1:nSess
    fprintf('  %-5s -> saved row %d\n', ...
        mouseIDs{s},savedRowForRequestedAnimal(s));
end

%% ---------------- behavior names ----------------

canonicalNames = { ...
    'climbdown', ...
    'climbup', ...
    'eating', ...
    'grooming', ...
    'jumpdown', ...
    'jumping', ...
    'rearing', ...
    'still', ...
    'walkflat', ...
    'walkgrid'};

behaviorLabels = strings(1,nBeh);

for bIdx = 1:nBeh
    beh = behaviors(bIdx);

    if beh >= 1 && beh <= 10
        behaviorLabels(bIdx) = sprintf('%d: %s',beh,canonicalNames{beh});
    else
        behaviorLabels(bIdx) = sprintf('Behavior %d',beh);
    end
end

%% ========================================================================
%% collect values ONCE
%% ========================================================================

peakLagMat = nan(nSess,nBeh);
peakCorrMat = nan(nSess,nBeh);

lagCILowMat = nan(nSess,nBeh);
lagCIHighMat = nan(nSess,nBeh);

nTrialsMat = nan(nSess,nBeh);
shiftThreshMat = nan(nSess,nBeh);
nAcceptedPermMat = zeros(nSess,nBeh);

permLagCell = cell(nSess,nBeh);

for bIdx = 1:nBeh

    B = behaviorResults{bIdx};

    if isempty(B)
        continue;
    end

    for s = 1:nSess

        savedRow = savedRowForRequestedAnimal(s);

        if isfield(B,'real_peakLagSec') && numel(B.real_peakLagSec) >= savedRow
            peakLagMat(s,bIdx) = double(B.real_peakLagSec(savedRow));
        end

        if isfield(B,'real_peakCorr') && numel(B.real_peakCorr) >= savedRow
            peakCorrMat(s,bIdx) = double(B.real_peakCorr(savedRow));
        end

        if isfield(B,'real_nTrialsUsed') && numel(B.real_nTrialsUsed) >= savedRow
            nTrialsMat(s,bIdx) = double(B.real_nTrialsUsed(savedRow));
        end

        if isfield(B,'shiftCorrUpper95') && numel(B.shiftCorrUpper95) >= savedRow
            shiftThreshMat(s,bIdx) = double(B.shiftCorrUpper95(savedRow));
        end

        if ~isfield(B,'perm_peakLagSec') || ...
                size(B.perm_peakLagSec,1) < savedRow

            permLagCell{s,bIdx} = [];
            continue;
        end

        allPermLags = double(B.perm_peakLagSec(savedRow,:));

        % Safe accepted-permutation handling.
        if isfield(B,'perm_accepted') && ...
                size(B.perm_accepted,1) >= savedRow

            acceptedRaw = double(B.perm_accepted(savedRow,:));

            % Only an explicit 1 counts as accepted.
            acceptedMask = isfinite(acceptedRaw) & (acceptedRaw == 1);

            % Protect against dimension mismatch.
            nCommon = min(numel(acceptedMask),numel(allPermLags));
            acceptedMask = acceptedMask(1:nCommon);
            allPermLags = allPermLags(1:nCommon);

        else
            acceptedMask = true(size(allPermLags));
        end

        useMask = acceptedMask & isfinite(allPermLags);

        permLags = allPermLags(useMask);
        permLagCell{s,bIdx} = permLags(:)';

        nAcceptedPermMat(s,bIdx) = numel(permLags);

        % THIS is the single source of truth for both plots and tables.
        if numel(permLags) >= 2
            ci = prctile(permLags,[2.5 97.5]);
            lagCILowMat(s,bIdx) = ci(1);
            lagCIHighMat(s,bIdx) = ci(2);
        end
    end
end

%% ========================================================================
%% 1. XCORR CURVES
%% ========================================================================

for bIdx = 1:nBeh

    B = behaviorResults{bIdx};

    if isempty(B)
        continue;
    end

    figure( ...
        'Name',sprintf('Chunked popavg xcorr - %s',behaviorLabels(bIdx)), ...
        'Color','w');

    tiledlayout(2,ceil(nSess/2), ...
        'TileSpacing','compact', ...
        'Padding','compact');

    for s = 1:nSess

        nexttile;
        hold on;

        savedRow = savedRowForRequestedAnimal(s);

        if ~isfield(B,'real_xc') || isempty(B.real_xc) || ...
                size(B.real_xc,1) < savedRow || ...
                all(~isfinite(B.real_xc(savedRow,:))) || ...
                ~isfield(B,'lags') || isempty(B.lags)

            title(sprintf('%s - missing',animalLabels(s)));
            axis off;
            continue;
        end

        if isfield(B,'binSize') && isscalar(B.binSize) && isfinite(B.binSize)
            lagsSec = double(B.lags) * double(B.binSize);
        else
            lagsSec = double(B.lags);
        end

        plot(lagsSec,double(B.real_xc(savedRow,:)),'k','LineWidth',2);

        if isfinite(peakLagMat(s,bIdx))
            xline(peakLagMat(s,bIdx),'r-','LineWidth',1.5);
        end

        if isfinite(shiftThreshMat(s,bIdx))
            yline(shiftThreshMat(s,bIdx),'--', ...
                'Color',[0.2 0.2 0.2], ...
                'LineWidth',1.25);
        end

        xlabel('Lag (s)');
        ylabel('Correlation');

        title(sprintf('%s | n=%g trials', ...
            animalLabels(s),nTrialsMat(s,bIdx)));

        box off;
    end

    sgtitle(sprintf('Chunked trial-averaged cortex popavg: %s', ...
        behaviorLabels(bIdx)));
end

%% ========================================================================
%% 2. CORRECTED PEAK-LAG DOT PLOT + EXACT PERMUTATION 95% BOUNDS
%% ========================================================================

for bIdx = 1:nBeh

    figure( ...
        'Name',sprintf('Chunked peak lag summary - %s',behaviorLabels(bIdx)), ...
        'Color','w');

    hold on;

    xPos = 1:nSess;

    for s = 1:nSess

        lo = lagCILowMat(s,bIdx);
        hi = lagCIHighMat(s,bIdx);
        actual = peakLagMat(s,bIdx);

        if isfinite(lo) && isfinite(hi)

            % Main vertical 95% permutation interval.
            line([xPos(s) xPos(s)],[lo hi], ...
                'Color',[0.45 0.45 0.45], ...
                'LineWidth',2);

            % Explicit caps at EXACT table values.
            capHalfWidth = 0.09;

            line([xPos(s)-capHalfWidth xPos(s)+capHalfWidth],[lo lo], ...
                'Color',[0.45 0.45 0.45], ...
                'LineWidth',2);

            line([xPos(s)-capHalfWidth xPos(s)+capHalfWidth],[hi hi], ...
                'Color',[0.45 0.45 0.45], ...
                'LineWidth',2);

            % Small endpoint markers make exact percentile locations obvious.
            scatter(xPos(s),lo,28,[0.45 0.45 0.45],'filled','Marker','_');
            scatter(xPos(s),hi,28,[0.45 0.45 0.45],'filled','Marker','_');
        end

        % Plot the actual lag separately. It does NOT need to lie inside
        % the permutation interval.
        if isfinite(actual)
            scatter(xPos(s),actual,70,'k','filled');
        end
    end

    yline(0,'k:');

    xlim([0.5 nSess+0.5]);

    xticks(xPos);

    % Include the number of real behavior events/trials for each animal.
    trialTickLabels = strings(1,nSess);
    for s = 1:nSess
        if isfinite(nTrialsMat(s,bIdx))
            trialTickLabels(s) = sprintf('%s\nn=%g trials', ...
                animalLabels(s),nTrialsMat(s,bIdx));
        else
            trialTickLabels(s) = sprintf('%s\nn=NA',animalLabels(s));
        end
    end
    xticklabels(cellstr(trialTickLabels));

    xlabel('Animal');
    ylabel('Peak lag (s)');

    title(sprintf('%s: actual peak lag + exact 2.5th-97.5th permutation bounds', ...
        behaviorLabels(bIdx)));

    box off;
    grid on;
end

%% ========================================================================
%% 3. PER-ANIMAL PERMUTATION HISTOGRAMS
%% ========================================================================

for bIdx = 1:nBeh

    [commonEdges,commonXLim] = getCommonHistogramEdges( ...
        permLagCell(:,bIdx),peakLagMat(:,bIdx),24);

    figure( ...
        'Name',sprintf('Chunked permutation peak lags - %s',behaviorLabels(bIdx)), ...
        'Color','w');

    tiledlayout(2,ceil(nSess/2), ...
        'TileSpacing','compact', ...
        'Padding','compact');

    for s = 1:nSess

        nexttile;
        hold on;

        permLags = permLagCell{s,bIdx};

        if isempty(permLags)
            title(sprintf('%s - no accepted perms',animalLabels(s)));
            axis off;
            continue;
        end

        histogram(permLags, ...
            'BinEdges',commonEdges, ...
            'FaceColor',[0.3 0.6 0.8], ...
            'EdgeColor','none');

        % IMPORTANT: use exact stored CI values, do not recalculate.
        lo = lagCILowMat(s,bIdx);
        hi = lagCIHighMat(s,bIdx);

        if isfinite(lo)
            xline(lo,'--','Color',[0.2 0.2 0.2],'LineWidth',1.5);
        end

        if isfinite(hi)
            xline(hi,'--','Color',[0.2 0.2 0.2],'LineWidth',1.5);
        end

        if isfinite(peakLagMat(s,bIdx))
            xline(peakLagMat(s,bIdx),'r-','LineWidth',1.5);
        end

        xlim(commonXLim);

        xlabel('Peak lag (s)');
        ylabel('Count');

        if isfinite(nTrialsMat(s,bIdx))
            title(sprintf('%s | n=%g trials | %d accepted perms', ...
                animalLabels(s),nTrialsMat(s,bIdx),nAcceptedPermMat(s,bIdx)));
        else
            title(sprintf('%s | n=NA trials | %d accepted perms', ...
                animalLabels(s),nAcceptedPermMat(s,bIdx)));
        end

        box off;
    end

    sgtitle(sprintf('%s: accepted permutation peak-lag distributions', ...
        behaviorLabels(bIdx)));
end

%% ========================================================================
%% 4. COMBINED PERMUTATION HISTOGRAM
%% ========================================================================

for bIdx = 1:nBeh

    allPermLags = [];

    for s = 1:nSess
        allPermLags = [allPermLags,permLagCell{s,bIdx}]; %#ok<AGROW>
    end

    allPermLags = allPermLags(isfinite(allPermLags));

    if isempty(allPermLags)
        continue;
    end

    [commonEdges,commonXLim] = getCommonHistogramEdges( ...
        permLagCell(:,bIdx),peakLagMat(:,bIdx),24);

    figure( ...
        'Name',sprintf('Chunked combined permutations - %s',behaviorLabels(bIdx)), ...
        'Color','w');

    hold on;

    hHist = histogram(allPermLags, ...
        'BinEdges',commonEdges, ...
        'FaceColor',[0.3 0.6 0.8], ...
        'EdgeColor','none');

    prcAll = prctile(allPermLags,[2.5 97.5]);

    hLow = xline(prcAll(1),'--', ...
        'Color',[0.2 0.2 0.2], ...
        'LineWidth',1.5);

    hHigh = xline(prcAll(2),'--', ...
        'Color',[0.2 0.2 0.2], ...
        'LineWidth',1.5);

    co = lines(nSess);
    actualLines = gobjects(0,1);
    actualLabels = {};

    for s = 1:nSess
        if isfinite(peakLagMat(s,bIdx))
            actualLines(end+1,1) = xline( ...
                peakLagMat(s,bIdx),'-', ...
                'Color',co(s,:), ...
                'LineWidth',1.5); %#ok<AGROW>

            actualLabels{end+1} = sprintf( ...
                '%s actual = %.3g s', ...
                animalLabels(s),peakLagMat(s,bIdx)); %#ok<AGROW>
        end
    end

    xlim(commonXLim);

    xlabel('Peak lag (s)');
    ylabel('Count');

    title(sprintf('%s: combined permutation distribution', ...
        behaviorLabels(bIdx)));

    box off;

    legHandles = [hHist,hLow,hHigh,actualLines(:)'];

    legEntries = { ...
        'Accepted permuted lags', ...
        sprintf('2.5%% all = %.3g s',prcAll(1)), ...
        sprintf('97.5%% all = %.3g s',prcAll(2))};

    legEntries = [legEntries,actualLabels];

    legend(legHandles,legEntries,'Location','best','Box','off');
end

%% ========================================================================
%% 5. OVERALL ANIMAL x BEHAVIOR PEAK-LAG SUMMARY
%% ========================================================================

figure( ...
    'Name','Chunked popavg peak lag by animal and behavior', ...
    'Color','w');

imagesc(peakLagMat);
colorbar;

xticks(1:nBeh);
xticklabels(cellstr(behaviorLabels));
xtickangle(45);

yticks(1:nSess);
yticklabels(cellstr(animalLabels));

xlabel('Behavior');
ylabel('Animal');

title('Actual chunked trial-averaged population peak lag (s)');

%% ========================================================================
%% PRINT SUMMARY: EXACT SAME CI VALUES USED BY DOT PLOT
%% ========================================================================

fprintf('\n============================================================\n');
fprintf('CHUNKED TRIAL-AVERAGED POPAVG PER-BEHAVIOR PLOT SUMMARY\n');
fprintf('============================================================\n');

for bIdx = 1:nBeh

    fprintf('\n%s\n',behaviorLabels(bIdx));

    T = table( ...
        animalLabels(:), ...
        nTrialsMat(:,bIdx), ...
        peakLagMat(:,bIdx), ...
        peakCorrMat(:,bIdx), ...
        lagCILowMat(:,bIdx), ...
        lagCIHighMat(:,bIdx), ...
        nAcceptedPermMat(:,bIdx), ...
        shiftThreshMat(:,bIdx), ...
        'VariableNames',{ ...
            'Animal', ...
            'NTrials', ...
            'PeakLag_s', ...
            'PeakCorr', ...
            'PermCI_Low_s', ...
            'PermCI_High_s', ...
            'NAcceptedPerms', ...
            'ShiftCorrUpper95'});

    disp(T);
end

end

%% ========================================================================
%% helper: common histogram edges
%% ========================================================================

function [edges,xLim] = getCommonHistogramEdges( ...
    permCells,actualVals,nBins)

allPerm = [];

for i = 1:numel(permCells)
    if ~isempty(permCells{i})
        allPerm = [allPerm,permCells{i}(:)']; %#ok<AGROW>
    end
end

allPerm = allPerm(isfinite(allPerm));
actualVals = actualVals(isfinite(actualVals));

if isempty(allPerm) && isempty(actualVals)

    xLim = [-0.5 0.5];

else

    allVals = [allPerm(:);actualVals(:)];

    xMin = min(allVals);
    xMax = max(allVals);

    if xMin == xMax
        pad = 0.001;
    else
        pad = 0.05*(xMax-xMin);
    end

    xLim = [xMin-pad,xMax+pad];
end

edges = linspace(xLim(1),xLim(2),nBins+1);

end
