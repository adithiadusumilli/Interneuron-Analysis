function runChunkedPopavgPerBehavior_plotOnly_getMouseDataNames(mouseIDs, baseSessionNames, probeRegions, saveFile)
% Plot-only companion for the combined CHUNKED trial-averaged
% population-average behavior-specific analysis.

% Combined structure:
%   allSessions.behaviorResults{bIdx}

% Each behavior result contains all 6 animals at once:
%   B.real_xc                 [nAnimals x nLags]
%   B.real_peakLagSec         [nAnimals x 1]
%   B.real_peakCorr           [nAnimals x 1]
%   B.real_nTrialsUsed        [nAnimals x 1]
%   B.perm_peakLagSec         [nAnimals x 300]
%   B.perm_accepted           [nAnimals x 300]
%   B.shiftCorrUpper95        [nAnimals x 1]

% Plot logic parallels the non-chunked version:
%   1) For each behavior: xcorr curves for all animals
%   2) For each behavior: actual peak lag + 95% permutation CI across animals
%   3) For each behavior: per-animal permutation peak-lag histograms
%   4) For each behavior: combined permutation histogram across animals
%   5) Heatmap summary of actual peak lags (animal x behavior)

% Chunked-specific differences:
%   - permutation CIs/histograms use ACCEPTED permutations only
%   - the animal x behavior-specific shiftCorrUpper95 threshold is drawn
%     on each xcorr panel as a horizontal dashed line
%   - trial counts are shown instead of continuous-timepoint counts

% Default combined file location outside Quest:
%   C:\Users\mirilab\Documents\GlobusTransfer\concatCrossCorrPerCanonicalBehavior_classifier_chunked_trialavg_popavg_ALL_SESSIONS.mat

% Example:
%   runChunkedPopavgPerBehavior_plotOnly_getMouseDataNames

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
    error('Saved combined chunked results file not found:\n%s', saveFile);
end

S = load(saveFile,'allSessions');

if ~isfield(S,'allSessions')
    error('File does not contain variable "allSessions":\n%s',saveFile);
end

allSessions = S.allSessions;

if ~isfield(allSessions,'behaviorResults') || isempty(allSessions.behaviorResults)
    error(['allSessions.behaviorResults is missing or empty. ' ...
           'This plotter expects the v2 combined structure.']);
end

behaviorResults = allSessions.behaviorResults;

if isfield(allSessions,'behaviors') && ~isempty(allSessions.behaviors)
    behaviors = allSessions.behaviors;
else
    behaviors = nan(1,numel(behaviorResults));
    for bIdx = 1:numel(behaviorResults)
        behaviors(bIdx) = behaviorResults{bIdx}.beh;
    end
end

nBeh = numel(behaviors);

%% ---------------- animal labels / sanity check ----------------

animalLabels = strings(1,nSess);

for s = 1:nSess
    animalLabels(s) = string(mouseIDs{s});

    if isfield(allSessions,'mouseIDs') && numel(allSessions.mouseIDs) >= s
        savedMouse = string(allSessions.mouseIDs{s});

        if savedMouse ~= string(mouseIDs{s})
            warning(['Animal %d saved mouseID is %s, but supplied mouseID is %s. ' ...
                     'Plots use supplied ordering.'], ...
                     s,savedMouse,string(mouseIDs{s}));
        end
    end

    if ~isstruct(resolvedDataNames{s}) || ...
            ~isfield(resolvedDataNames{s},'processedDataFolder') || ...
            isempty(resolvedDataNames{s}.processedDataFolder)

        warning('getMouseDataNames returned no processedDataFolder for %s.', ...
            mouseIDs{s});
    end
end

%% ---------------- canonical behavior names ----------------

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
        behaviorLabels(bIdx) = sprintf('%d: %s', ...
            beh,canonicalNames{beh});
    else
        behaviorLabels(bIdx) = sprintf('Behavior %d',beh);
    end
end

%% ---------------- collect matrices / cells ----------------

peakLagMat = nan(nSess,nBeh);
peakCorrMat = nan(nSess,nBeh);

lagCILowMat = nan(nSess,nBeh);
lagCIHighMat = nan(nSess,nBeh);

nTrialsMat = nan(nSess,nBeh);
shiftThreshMat = nan(nSess,nBeh);
nAcceptedPermMat = nan(nSess,nBeh);

permLagCell = cell(nSess,nBeh);

for bIdx = 1:nBeh

    B = behaviorResults{bIdx};

    if isempty(B)
        continue;
    end

    for s = 1:nSess

        if isfield(B,'real_peakLagSec') && numel(B.real_peakLagSec) >= s
            peakLagMat(s,bIdx) = B.real_peakLagSec(s);
        end

        if isfield(B,'real_peakCorr') && numel(B.real_peakCorr) >= s
            peakCorrMat(s,bIdx) = B.real_peakCorr(s);
        end

        if isfield(B,'real_nTrialsUsed') && numel(B.real_nTrialsUsed) >= s
            nTrialsMat(s,bIdx) = B.real_nTrialsUsed(s);
        end

        if isfield(B,'shiftCorrUpper95') && numel(B.shiftCorrUpper95) >= s
            shiftThreshMat(s,bIdx) = B.shiftCorrUpper95(s);
        end

        if isfield(B,'perm_peakLagSec') && size(B.perm_peakLagSec,1) >= s

            allPermLags = B.perm_peakLagSec(s,:);

            if isfield(B,'perm_accepted') && size(B.perm_accepted,1) >= s
                acceptedMask = logical(B.perm_accepted(s,:));
            else
                acceptedMask = isfinite(allPermLags);
            end

            useMask = acceptedMask & isfinite(allPermLags);

            permLags = allPermLags(useMask);

            permLagCell{s,bIdx} = permLags(:)';
            nAcceptedPermMat(s,bIdx) = numel(permLags);

            if numel(permLags) >= 2
                ci = prctile(permLags,[2.5 97.5]);

                lagCILowMat(s,bIdx) = ci(1);
                lagCIHighMat(s,bIdx) = ci(2);
            end
        else
            permLagCell{s,bIdx} = [];
            nAcceptedPermMat(s,bIdx) = 0;
        end
    end
end

%% ========================================================================
%% 1. XCORR CURVES: one figure per behavior, all animals
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

        if ~isfield(B,'real_xc') || ...
                isempty(B.real_xc) || ...
                size(B.real_xc,1) < s || ...
                all(~isfinite(B.real_xc(s,:))) || ...
                ~isfield(B,'lags') || ...
                isempty(B.lags)

            title(sprintf('%s - missing',animalLabels(s)));
            axis off;
            continue;
        end

        binSize = NaN;

        if isfield(B,'binSize') && isscalar(B.binSize)
            binSize = B.binSize;
        end

        if isfinite(binSize)
            lagsSec = B.lags * binSize;
        else
            lagsSec = B.lags;
        end

        plot(lagsSec,B.real_xc(s,:),'k','LineWidth',2);

        if isfinite(peakLagMat(s,bIdx))
            xline(peakLagMat(s,bIdx),'r-','LineWidth',1.5);
        end

        % Behavior-specific permutation correlation threshold.
        if isfinite(shiftThreshMat(s,bIdx))
            yline(shiftThreshMat(s,bIdx),'--', ...
                'Color',[0.2 0.2 0.2], ...
                'LineWidth',1.25);
        end

        xlabel('Lag (s)');
        ylabel('Correlation');

        title(sprintf('%s | n=%g trials', ...
            animalLabels(s), ...
            nTrialsMat(s,bIdx)));

        box off;
    end

    sgtitle(sprintf('Chunked trial-averaged cortex popavg: %s', ...
        behaviorLabels(bIdx)));
end

%% ========================================================================
%% 2. PEAK-LAG SUMMARY + ACCEPTED-PERMUTATION CI
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

        if isfinite(lo) && isfinite(hi)
            line([xPos(s) xPos(s)],[lo hi], ...
                'Color',[0.6 0.6 0.6], ...
                'LineWidth',2);
        end
    end

    scatter(xPos,peakLagMat(:,bIdx),70,'k','filled');
    yline(0,'k:');

    xlim([0.5 nSess+0.5]);

    xticks(xPos);
    xticklabels(cellstr(animalLabels));

    xlabel('Animal');
    ylabel('Peak lag (s)');

    title(sprintf('%s: actual peak lag with 95%% permutation CI', ...
        behaviorLabels(bIdx)));

    box off;
    grid on;
end

%% ========================================================================
%% 3. PER-ANIMAL ACCEPTED-PERMUTATION HISTOGRAMS
%% ========================================================================

for bIdx = 1:nBeh

    [commonEdges,commonXLim] = ...
        getCommonHistogramEdges( ...
            permLagCell(:,bIdx), ...
            peakLagMat(:,bIdx), ...
            24);

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

        prcLag = prctile(permLags,[2.5 97.5]);

        xline(prcLag(1),'--', ...
            'Color',[0.2 0.2 0.2], ...
            'LineWidth',1.5);

        xline(prcLag(2),'--', ...
            'Color',[0.2 0.2 0.2], ...
            'LineWidth',1.5);

        if isfinite(peakLagMat(s,bIdx))
            xline(peakLagMat(s,bIdx),'r-','LineWidth',1.5);
        end

        xlim(commonXLim);

        xlabel('Peak lag (s)');
        ylabel('Count');

        title(sprintf('%s | %d accepted perms', ...
            animalLabels(s), ...
            numel(permLags)));

        box off;
    end

    sgtitle(sprintf('%s: accepted permutation peak-lag distributions', ...
        behaviorLabels(bIdx)));
end

%% ========================================================================
%% 4. COMBINED ACCEPTED-PERMUTATION HISTOGRAM
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

    [commonEdges,commonXLim] = ...
        getCommonHistogramEdges( ...
            permLagCell(:,bIdx), ...
            peakLagMat(:,bIdx), ...
            24);

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
                peakLagMat(s,bIdx), ...
                '-', ...
                'Color',co(s,:), ...
                'LineWidth',1.5); %#ok<AGROW>

            actualLabels{end+1} = sprintf( ...
                '%s actual = %.3g s', ...
                animalLabels(s), ...
                peakLagMat(s,bIdx)); %#ok<AGROW>
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

    legend( ...
        legHandles, ...
        legEntries, ...
        'Location','best', ...
        'Box','off');
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
%% PRINT SUMMARY
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
    permCells, ...
    actualVals, ...
    nBins)

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
