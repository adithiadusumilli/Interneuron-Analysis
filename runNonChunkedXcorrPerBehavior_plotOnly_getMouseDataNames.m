function runNonChunkedXcorrPerBehavior_plotOnly_getMouseDataNames(mouseIDs, baseSessionNames, probeRegions, saveFile)
% Plot-only for: runNonChunkedXcorrPerBehavior_saveOnly_getMouseDataNames

% Uses the saved non-chunked, population-averaged, behavior-specific results.
% No trial averaging is performed here; this function only regenerates plots.

% Inputs:
%   mouseIDs - cell array of mouse IDs
%   baseSessionNames - cell array of base session names
%   probeRegions - cell array of probe regions
%   saveFile - optional path to saved results .mat
%
% If saveFile is omitted, uses:
%   X:\David\AnalysesData\nonchunked_xcorr_by_classifier_cortex_all6Sessions_saved.mat

% For consistency with the save-only analysis, getMouseDataNames is called
% for each supplied animal/session/probe combination and checked against the
% metadata saved in results.

% Plots:
%   1) For each behavior: xcorr curves for all animals
%   2) For each behavior: actual peak lag + 95% permutation CI across animals
%   3) For each behavior: per-animal permutation peak-lag histograms
%   4) For each behavior: combined permutation histogram across animals
%   5) Heatmap-style summary of actual peak lags (animal x behavior)

% Positive/negative lag interpretation is inherited directly from the
% save-only function's computeManualXCorrVec convention.

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
        'X:\David\AnalysesData', ...
        'nonchunked_xcorr_by_classifier_cortex_all6Sessions_saved.mat');
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

%% ---------------- load saved analysis ----------------

if ~isfile(saveFile)
    error('Saved results file not found:\n%s', saveFile);
end

S = load(saveFile, 'results');

if ~isfield(S,'results')
    error('File does not contain variable "results":\n%s', saveFile);
end

results = S.results;

if ~isfield(results,'sessions') || isempty(results.sessions)
    error('results.sessions is missing or empty.');
end

if isfield(results,'behaviors') && ~isempty(results.behaviors)
    behaviors = results.behaviors;
else
    % Recover behavior values from the first nonempty session.
    behaviors = [];
    for s = 1:numel(results.sessions)
        if isfield(results.sessions(s),'beh') && ~isempty(results.sessions(s).beh)
            behaviors = [results.sessions(s).beh.beh];
            break;
        end
    end
end

if isempty(behaviors)
    error('Could not determine behavior values from saved results.');
end

nBeh = numel(behaviors);

%% ---------------- animal labels ----------------

animalLabels = strings(1,nSess);

for s = 1:nSess
    animalLabels(s) = string(mouseIDs{s});

    % Sanity check saved ordering against requested ordering.
    if s <= numel(results.sessions) && ...
            isfield(results.sessions(s),'mouseID') && ...
            ~isempty(results.sessions(s).mouseID)

        savedMouse = string(results.sessions(s).mouseID);

        if savedMouse ~= string(mouseIDs{s})
            warning(['Session %d saved mouseID is %s, but supplied mouseID is %s. ' ...
                     'Plots use the supplied mouse ordering.'], ...
                     s, savedMouse, string(mouseIDs{s}));
        end
    end

    % Sanity check processed folder when available.
    if s <= numel(results.sessions) && ...
            isfield(results.sessions(s),'processedDataFolder') && ...
            ~isempty(results.sessions(s).processedDataFolder) && ...
            isfield(resolvedDataNames{s},'processedDataFolder') && ...
            ~isempty(resolvedDataNames{s}.processedDataFolder)

        savedDir = string(results.sessions(s).processedDataFolder);
        resolvedDir = string(resolvedDataNames{s}.processedDataFolder);

        if ~strcmpi(savedDir,resolvedDir)
            warning(['%s: saved processedDataFolder differs from current ' ...
                     'getMouseDataNames output.\nSaved: %s\nCurrent: %s'], ...
                     mouseIDs{s}, savedDir, resolvedDir);
        end
    end
end

%% ---------------- behavior names ----------------

canonicalNames = { ...
    'unlabeled', ...
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

    if beh >= 0 && beh <= 10
        behaviorLabels(bIdx) = sprintf('%d: %s', ...
            beh, canonicalNames{beh+1});
    else
        behaviorLabels(bIdx) = sprintf('Behavior %d', beh);
    end
end

%% ---------------- collect matrices/cells ----------------

peakLagMat = nan(nSess,nBeh);
peakCorrMat = nan(nSess,nBeh);
lagCILowMat = nan(nSess,nBeh);
lagCIHighMat = nan(nSess,nBeh);
nTimepointsMat = nan(nSess,nBeh);
permLagCell = cell(nSess,nBeh);

for s = 1:nSess

    if s > numel(results.sessions) || ...
            ~isfield(results.sessions(s),'beh') || ...
            isempty(results.sessions(s).beh)
        continue;
    end

    sessBeh = results.sessions(s).beh;

    for bIdx = 1:nBeh
        beh = behaviors(bIdx);

        thisIdx = find([sessBeh.beh] == beh,1,'first');

        if isempty(thisIdx)
            continue;
        end

        B = sessBeh(thisIdx);

        peakLagMat(s,bIdx) = scalarOrNaN(B,'peakLagSec');
        peakCorrMat(s,bIdx) = scalarOrNaN(B,'peakCorr');
        nTimepointsMat(s,bIdx) = scalarOrNaN(B,'nTimepoints');

        if isfield(B,'lagCI') && numel(B.lagCI) >= 2
            lagCILowMat(s,bIdx) = B.lagCI(1);
            lagCIHighMat(s,bIdx) = B.lagCI(2);
        end

        if isfield(B,'permPeakLags') && ~isempty(B.permPeakLags)
            x = B.permPeakLags(:)';
            permLagCell{s,bIdx} = x(isfinite(x));
        else
            permLagCell{s,bIdx} = [];
        end
    end
end

%% ========================================================================
%% 1. XCORR CURVES: one figure per behavior, all animals
%% ========================================================================

for bIdx = 1:nBeh

    beh = behaviors(bIdx);

    figure( ...
        'Name',sprintf('Non-chunked popavg xcorr - %s',behaviorLabels(bIdx)), ...
        'Color','w');

    tiledlayout(2,ceil(nSess/2), ...
        'TileSpacing','compact', ...
        'Padding','compact');

    for s = 1:nSess

        nexttile;
        hold on;

        B = getBehaviorStruct(results,s,beh);

        if isempty(B) || ...
                ~isfield(B,'xc') || isempty(B.xc) || ...
                ~isfield(B,'lagsSec') || isempty(B.lagsSec)

            title(sprintf('%s - missing',animalLabels(s)));
            axis off;
            continue;
        end

        plot(B.lagsSec,B.xc,'k','LineWidth',2);

        if isfield(B,'peakLagSec') && isfinite(B.peakLagSec)
            xline(B.peakLagSec,'r-','LineWidth',1.5);
        end

        % ctrlCorrCI is a correlation threshold/interval, so it belongs
        % on the y-axis of the xcorr plot.
        if isfield(B,'ctrlCorrCI') && numel(B.ctrlCorrCI) >= 2
            if isfinite(B.ctrlCorrCI(1))
                yline(B.ctrlCorrCI(1),'--', ...
                    'Color',[0.2 0.2 0.2], ...
                    'LineWidth',1);
            end
            if isfinite(B.ctrlCorrCI(2))
                yline(B.ctrlCorrCI(2),'--', ...
                    'Color',[0.2 0.2 0.2], ...
                    'LineWidth',1);
            end
        end

        xlabel('Lag (s)');
        ylabel('Correlation');
        title(sprintf('%s | n=%g timepoints', ...
            animalLabels(s),nTimepointsMat(s,bIdx)));

        box off;
    end

    sgtitle(sprintf('Non-chunked cortex popavg: %s',behaviorLabels(bIdx)));
end

%% ========================================================================
%% 2. PEAK-LAG SUMMARY + PERMUTATION CI: one figure per behavior
%% ========================================================================

for bIdx = 1:nBeh

    figure( ...
        'Name',sprintf('Peak lag summary - %s',behaviorLabels(bIdx)), ...
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
%% 3. PER-ANIMAL PERMUTATION HISTOGRAMS: one figure per behavior
%% ========================================================================

for bIdx = 1:nBeh

    [commonEdges,commonXLim] = ...
        getCommonHistogramEdges( ...
            permLagCell(:,bIdx), ...
            peakLagMat(:,bIdx), ...
            24);

    figure( ...
        'Name',sprintf('Permutation peak lags - %s',behaviorLabels(bIdx)), ...
        'Color','w');

    tiledlayout(2,ceil(nSess/2), ...
        'TileSpacing','compact', ...
        'Padding','compact');

    for s = 1:nSess

        nexttile;
        hold on;

        permLags = permLagCell{s,bIdx};

        if isempty(permLags)
            title(sprintf('%s - no perms',animalLabels(s)));
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
        title(sprintf('%s | %d perms',animalLabels(s),numel(permLags)));

        box off;
    end

    sgtitle(sprintf('%s: permutation peak-lag distributions', ...
        behaviorLabels(bIdx)));
end

%% ========================================================================
%% 4. COMBINED PERMUTATION HISTOGRAM: one figure per behavior
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
        'Name',sprintf('Combined permutations - %s',behaviorLabels(bIdx)), ...
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
        'Permuted lags', ...
        sprintf('2.5%% all = %.3g s',prcAll(1)), ...
        sprintf('97.5%% all = %.3g s',prcAll(2))};

    legEntries = [legEntries,actualLabels];

    legend(legHandles,legEntries, ...
        'Location','best', ...
        'Box','off');
end

%% ========================================================================
%% 5. OVERALL ANIMAL x BEHAVIOR PEAK-LAG SUMMARY
%% ========================================================================

figure( ...
    'Name','Non-chunked popavg peak lag by animal and behavior', ...
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
title('Actual non-chunked population peak lag (s)');

%% ========================================================================
%% PRINT SUMMARY
%% ========================================================================

fprintf('\n============================================================\n');
fprintf('NON-CHUNKED POPAVG PER-BEHAVIOR PLOT SUMMARY\n');
fprintf('============================================================\n');

for bIdx = 1:nBeh

    fprintf('\n%s\n',behaviorLabels(bIdx));

    T = table( ...
        animalLabels(:), ...
        nTimepointsMat(:,bIdx), ...
        peakLagMat(:,bIdx), ...
        peakCorrMat(:,bIdx), ...
        lagCILowMat(:,bIdx), ...
        lagCIHighMat(:,bIdx), ...
        'VariableNames',{ ...
            'Animal', ...
            'NTimepoints', ...
            'PeakLag_s', ...
            'PeakCorr', ...
            'PermCI_Low_s', ...
            'PermCI_High_s'});

    disp(T);
end

end

%% ========================================================================
%% helper: retrieve behavior struct
%% ========================================================================

function B = getBehaviorStruct(results,sessIdx,beh)

B = [];

if sessIdx > numel(results.sessions)
    return;
end

if ~isfield(results.sessions(sessIdx),'beh') || ...
        isempty(results.sessions(sessIdx).beh)
    return;
end

sessBeh = results.sessions(sessIdx).beh;
idx = find([sessBeh.beh] == beh,1,'first');

if isempty(idx)
    return;
end

B = sessBeh(idx);

end

%% ========================================================================
%% helper: safe scalar field
%% ========================================================================

function x = scalarOrNaN(S,fieldName)

x = NaN;

if isfield(S,fieldName) && ~isempty(S.(fieldName))
    val = S.(fieldName);

    if isscalar(val)
        x = val;
    end
end

end

%% ========================================================================
%% helper: common histogram edges
%% ========================================================================

function [edges,xLim] = getCommonHistogramEdges(permCells,actualVals,nBins)

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
