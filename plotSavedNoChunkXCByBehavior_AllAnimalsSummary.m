function plotSavedNoChunkXCByBehavior_AllAnimalsSummary(noChunkFile)
% plots saved nonchunked popavg xc results by behavior

% outputs:
%   1. four separate 2x5 lag-vs-correlation figures
%        - one figure per animal
%        - each subpanel = one behavior

%   2. one summary 2x5 figure
%        - each subpanel = one behavior
%        - each panel contains all animals:
%            black dot = actual peak lag
%            slightly offset gray vertical line = null 95% range

%   3. four separate 2x5 permutation histogram figures
%        - one figure per animal
%        - each subpanel = one behavior

% this version:
%   - plots all 10 behaviors for nonchunked results as long as data exist
%   - uses one overarching title per tiled figure
%   - uses shared legend beneath each tiled figure
%   - uses common x-axis limits within each animal's lag-vs-correlation figure
%   - uses common x-axis limits and bin edges within each animal's permutation histogram figure
%   - uses animal names only as titles where appropriate

% j run:
% plotSavedNoChunkXCByBehavior_AllAnimalsSummary("X:\David\AnalysesData\nonchunked_xcorr_by_classifier_cortex_allSessions_saved.mat")

arguments
    noChunkFile (1,1) string
end

behNums = 1:10;
behNames = {'climbdown','climbup','eating','grooming','jumpdown','jumping','rearing','still','walkflat','walkgrid'};

origColor = [0 0 0];
corrCIColor = [0 0.2 0.6];
peakLagColor = [0.95 0.45 0.35];
lagCIColor = [0 0 0];
permHistColor = [0.3 0.6 0.8];

titleFont = 16;
panelTitleFont = 15;
labelFont = 14;
tickFont = 13;
legendFont = 12;

A = load(noChunkFile, 'results');
R = A.results;

nSess = numel(R.sessions);
animalIDs = cell(1, nSess);

for i = 1:nSess
    a = regexp(R.sessions(i).baseDir, 'D\d+', 'match', 'once');
    if isempty(a)
        a = sprintf('Session %d', i);
    end
    animalIDs{i} = a;
end

binSize = 0.001;

% summary storage: session x behavior
actualLagMat = nan(nSess, numel(behNums));
lagCIMatLo = nan(nSess, numel(behNums));
lagCIMatHi = nan(nSess, numel(behNums));
boutMat = nan(nSess, numel(behNums));
durMat = nan(nSess, numel(behNums));

fprintf('\n================ nonchunked behavior plots ================\n')

%% ---- precompute common x-limits for lag-vs-correlation per animal ----
lagXLimBySess = nan(nSess, 2);
for s = 1:nSess
    allLagVals = [];
    for k = 1:numel(behNums)
        b = behNums(k);
        thisBeh = R.sessions(s).beh(b);

        if isempty(thisBeh) || ~isfield(thisBeh, 'lagsSec') || isempty(thisBeh.lagsSec)
            continue
        end

        allLagVals = [allLagVals; thisBeh.lagsSec(:)]; %#ok<AGROW>
    end

    allLagVals = allLagVals(~isnan(allLagVals));
    if ~isempty(allLagVals)
        lagXLimBySess(s,:) = [min(allLagVals) max(allLagVals)];
    else
        lagXLimBySess(s,:) = [-0.5 0.5];
    end
end

%% ---- precompute common x-limits and bin edges for permutation histograms per animal ----
permXLimBySess = nan(nSess, 2);
permEdgesBySess = cell(nSess,1);
nHistBins = 24;

for s = 1:nSess
    allPermVals = [];
    allActualVals = [];

    for k = 1:numel(behNums)
        b = behNums(k);
        thisBeh = R.sessions(s).beh(b);

        if isempty(thisBeh)
            continue
        end

        if isfield(thisBeh, 'permPeakLags') && ~isempty(thisBeh.permPeakLags)
            tmp = thisBeh.permPeakLags(:);
            tmp = tmp(~isnan(tmp));
            allPermVals = [allPermVals; tmp]; %#ok<AGROW>
        end

        if isfield(thisBeh, 'peakLagSec') && ~isempty(thisBeh.peakLagSec) && ~isnan(thisBeh.peakLagSec)
            allActualVals = [allActualVals; thisBeh.peakLagSec]; %#ok<AGROW>
        end
    end

    allVals = [allPermVals; allActualVals];
    allVals = allVals(~isnan(allVals));

    if isempty(allVals)
        xMin = -0.02;
        xMax = 0.02;
    else
        xMin = min(allVals);
        xMax = max(allVals);

        if xMin == xMax
            xPad = 0.001;
        else
            xPad = 0.05 * (xMax - xMin);
        end

        xMin = xMin - xPad;
        xMax = xMax + xPad;
    end

    permXLimBySess(s,:) = [xMin xMax];
    permEdgesBySess{s} = linspace(xMin, xMax, nHistBins + 1);
end

%% ---- plot 1: per animal 2x5 lag vs correlation ----
for s = 1:nSess
    figure('Name', sprintf('%s nonchunk lag vs corr by behavior', animalIDs{s}), 'Color', 'w');
    tile_lay1 = tiledlayout(2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tile_lay1, sprintf('%s', animalIDs{s}), 'FontSize', titleFont, 'FontWeight', 'bold');

    firstLegendDone = false;

    for k = 1:numel(behNums)
        b = behNums(k);
        thisBeh = R.sessions(s).beh(b);

        nexttile(tile_lay1, k); hold on;

        if isempty(thisBeh) || ~isfield(thisBeh, 'timeIdx') || isempty(thisBeh.timeIdx) || ...
                ~isfield(thisBeh, 'xc') || isempty(thisBeh.xc) || ...
                ~isfield(thisBeh, 'lagsSec') || isempty(thisBeh.lagsSec)
            title(behNames{k}, 'FontSize', panelTitleFont);
            axis off;
            continue
        end

        timeIdx = sort(thisBeh.timeIdx(:));
        if isempty(timeIdx)
            nBouts = 0;
        else
            nBouts = 1 + sum(diff(timeIdx) > 1);
        end
        durSec = numel(timeIdx) * binSize;

        lagsSec = thisBeh.lagsSec(:)';
        xc = thisBeh.xc(:)';
        peakLag = thisBeh.peakLagSec;
        corrCI = thisBeh.ctrlCorrCI;
        lagCI = thisBeh.lagCI;

        actualLagMat(s,k) = peakLag;
        if ~isempty(lagCI) && numel(lagCI) == 2 && ~any(isnan(lagCI))
            lagCIMatLo(s,k) = lagCI(1);
            lagCIMatHi(s,k) = lagCI(2);
        end
        boutMat(s,k) = nBouts;
        durMat(s,k) = durSec;

        hOrig = plot(lagsSec, xc, 'Color', origColor, 'LineWidth', 2);

        hCorr25 = gobjects(1);
        hCorr97 = gobjects(1);
        if ~isempty(corrCI) && ~any(isnan(corrCI))
            hCorr25 = yline(corrCI(1), '--', 'Color', corrCIColor, 'LineWidth', 1.4);
            hCorr97 = yline(corrCI(2), '--', 'Color', corrCIColor, 'LineWidth', 1.4);
        end

        hPeak = xline(peakLag, '--', 'Color', peakLagColor, 'LineWidth', 1.8);

        hLag1 = gobjects(1);
        if ~isempty(lagCI) && ~any(isnan(lagCI))
            hLag1 = xline(lagCI(1), '--', 'Color', lagCIColor, 'LineWidth', 1.4);
            xline(lagCI(2), '--', 'Color', lagCIColor, 'LineWidth', 1.4);
        end

        xlim(lagXLimBySess(s,:));
        xlabel('Lag (Seconds)', 'FontSize', labelFont);
        ylabel('Correlation', 'FontSize', labelFont);
        title(behNames{k}, 'FontSize', panelTitleFont);
        box off;
        set(gca, 'FontSize', tickFont, 'LineWidth', 1, 'TickDir', 'out');

        if ~firstLegendDone
            legHandles = [hOrig];
            legLabels = {'Real Cross-Correlation'};
            if isgraphics(hCorr25), legHandles(end+1) = hCorr25; legLabels{end+1} = '2.5% Shift Control Correlation'; end
            if isgraphics(hCorr97), legHandles(end+1) = hCorr97; legLabels{end+1} = '97.5% Shift Control Correlation'; end
            legHandles(end+1) = hPeak; legLabels{end+1} = 'Actual Peak Lag';
            if isgraphics(hLag1), legHandles(end+1) = hLag1; legLabels{end+1} = 'Lag 95% CI (Permutation)'; end

            lgd = legend(legHandles, legLabels, 'Orientation', 'horizontal');
            lgd.Layout.Tile = 'south';
            lgd.FontSize = legendFont;
            lgd.Box = 'off';
            firstLegendDone = true;
        end
    end
end

%% ---- plot 3: per animal 2x5 permutation histograms ----
for s = 1:nSess
    figure('Name', sprintf('%s nonchunk permutation histograms by behavior', animalIDs{s}), 'Color', 'w');
    tile_lay2 = tiledlayout(2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tile_lay2, sprintf('%s', animalIDs{s}), 'FontSize', titleFont, 'FontWeight', 'bold');

    sharedLegendHandles = gobjects(4,1);
    legendSet = false;

    for k = 1:numel(behNums)
        b = behNums(k);
        thisBeh = R.sessions(s).beh(b);

        nexttile(tile_lay2, k); hold on;

        if isempty(thisBeh) || ~isfield(thisBeh, 'timeIdx') || isempty(thisBeh.timeIdx) || ...
                ~isfield(thisBeh, 'permPeakLags') || isempty(thisBeh.permPeakLags)
            title(behNames{k}, 'FontSize', panelTitleFont);
            axis off;
            continue
        end

        permLags = thisBeh.permPeakLags(:);
        permLags = permLags(~isnan(permLags));

        if isempty(permLags)
            title(behNames{k}, 'FontSize', panelTitleFont);
            axis off;
            continue
        end

        hHist = histogram(permLags, ...
            'BinEdges', permEdgesBySess{s}, ...
            'FaceColor', permHistColor, ...
            'EdgeColor', 'none');

        prcLag = prctile(permLags, [2.5 97.5]);
        hCI1 = xline(prcLag(1), '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
        hCI2 = xline(prcLag(2), '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
        hPeak = xline(thisBeh.peakLagSec, 'r-', 'LineWidth', 1.8);

        if ~legendSet
            sharedLegendHandles(1) = hPeak;

            % dummy patch so blue appears correctly in legend
            hHistLegend = patch(nan, nan, permHistColor, 'EdgeColor', 'none');
            sharedLegendHandles(2) = hHistLegend;

            sharedLegendHandles(3) = hCI1;
            sharedLegendHandles(4) = hCI2;
            legendSet = true;
        end

        xlim(permXLimBySess(s,:));
        xlabel('Peak Lag (s)', 'FontSize', labelFont);
        ylabel('Count', 'FontSize', labelFont);
        title(behNames{k}, 'FontSize', panelTitleFont);
        box off;
        set(gca, 'FontSize', tickFont, 'LineWidth', 1, 'TickDir', 'out');
    end

    if legendSet
        lgd = legend(sharedLegendHandles, ...
            {'Actual Peak Lag', ...
             'Permuted Peak Lags', ...
             '2.5% Permutation Control Lag', ...
             '97.5% Permutation Control Lag'}, ...
            'Orientation', 'horizontal');
        lgd.Layout.Tile = 'south';
        lgd.FontSize = legendFont;
        lgd.Box = 'off';
    end
end

%% ---- plot 2: one summary 2x5 figure with all animals in each behavior panel ----
figure('Name', 'Nonchunked summary by behavior', 'Color', 'w');
tile_lay3 = tiledlayout(2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tile_lay3, 'Actual Peak Lags with Permutation Control Range by Behavior', ...
    'FontSize', titleFont, 'FontWeight', 'bold');

xBase = 1:nSess;
xOffsets = linspace(-0.18, 0.18, nSess);

% common y-limits across summary panels
summaryVals = [actualLagMat(:); lagCIMatLo(:); lagCIMatHi(:)];
summaryVals = summaryVals(~isnan(summaryVals));

if isempty(summaryVals)
    summaryYLim = [-0.02 0.02];
else
    yMin = min(summaryVals);
    yMax = max(summaryVals);
    if yMin == yMax
        yPad = 0.001;
    else
        yPad = 0.05 * (yMax - yMin);
    end
    summaryYLim = [yMin - yPad, yMax + yPad];
end

for k = 1:numel(behNums)
    nexttile(tile_lay3, k); hold on;

    for s = 1:nSess
        if ~isnan(lagCIMatLo(s,k)) && ~isnan(lagCIMatHi(s,k))
            xci = xBase(s) + xOffsets(s);
            line([xci xci], [lagCIMatLo(s,k) lagCIMatHi(s,k)], ...
                'Color', [0.6 0.6 0.6], 'LineWidth', 2);
        end

        if ~isnan(actualLagMat(s,k))
            plot(xBase(s), actualLagMat(s,k), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
        end
    end

    yline(0, 'k:');
    ylim(summaryYLim);
    xlim([0.5 nSess + 0.5]);
    xticks(1:nSess);
    xticklabels(animalIDs);
    ylabel('Peak Lag (s)', 'FontSize', labelFont);
    title(behNames{k}, 'FontSize', panelTitleFont);
    box off
    grid on
    set(gca, 'FontSize', tickFont, 'LineWidth', 1, 'TickDir', 'out');
end

hActual = plot(nan, nan, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
hNull = line([nan nan], [nan nan], 'Color', [0.6 0.6 0.6], 'LineWidth', 2);

lgd = legend([hActual hNull], ...
    {'Actual Peak Lag', 'Permutation Control Range'}, ...
    'Orientation', 'horizontal');
lgd.Layout.Tile = 'south';
lgd.FontSize = legendFont;
lgd.Box = 'off';

end
