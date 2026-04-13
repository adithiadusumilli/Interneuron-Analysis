function plotSavedChunkXC_SelectedBehaviors_AllAnimalsSummary(chunkFile, minTrials)
% plots saved chunked popavg xc results for selected classifier behaviors

% selected behaviors: climb up (2), eating (3), walk grid (10)

% outputs:
%   1. four separate 1x3 lag-vs-correlation figures
%        - one figure per animal
%        - each subpanel = one selected behavior

%   2. one summary 1x3 figure
%        - each subpanel = one selected behavior
%        - each panel contains all animals:
%            black dot = actual peak lag
%            slightly offset gray vertical line = permutation control range

%   3. four separate 1x3 permutation histogram figures
%        - one figure per animal
%        - each subpanel = one selected behavior

% this version matches the plotting style used in nonchunk func (plotSavedNoChunkXCByBehavior_AllAnimalsSummary)

% j run: plotSavedChunkXC_SelectedBehaviors_AllAnimalsSummary("C:\Users\mirilab\Documents\GlobusTransfer\combined_allAnimals_concatCrossCorrPerCanonicalBehavior_classifier.mat", 25)

arguments
    chunkFile (1,1) string
    minTrials (1,1) double = 25
end

behNums = [2 3 10];
behNames = {'climbup','eating','walkgrid'};

% prettier labels for plot titles only
behNamesPretty = {'Climb Up','Eating','Walk Grid'};

origColor = [0 0 0];
corrCIColor = [0 0.2 0.6];
peakLagColor = [0.95 0.45 0.35];
lagCIColor = [0 0 0];
permHistColor = [0.3 0.6 0.8];

titleFont = 17;
panelTitleFont = 16;
labelFont = 15;
tickFont = 14;
legendFont = 14;

C = load(chunkFile, ...
    'animalNames', 'baseDirs', 'behaviors', 'lags', 'binSize', ...
    'all_xc_real', 'all_peakLagSec_real', 'all_nTrials_real', ...
    'all_peakLagSec_perm', 'all_xcZeroLag_shift');

nSess = numel(C.animalNames);
animalIDs = cell(1, nSess);

for i = 1:nSess
    a = regexp(C.baseDirs{i}, 'D\d+', 'match', 'once');
    if isempty(a)
        a = sprintf('Session %d', i);
    end
    animalIDs{i} = a;
end

actualLagMat = nan(nSess, numel(behNums));
lagCIMatLo = nan(nSess, numel(behNums));
lagCIMatHi = nan(nSess, numel(behNums));
trialMat = nan(nSess, numel(behNums));

fprintf('\n================ chunked selected behavior plots ================\n')

%% ---- precompute common x-limits for lag-vs-correlation per animal ----
lagXLimBySess = nan(nSess, 2);
for s = 1:nSess
    allLagVals = [];

    for k = 1:numel(behNums)
        b = behNums(k);
        bIdx = find(C.behaviors == b, 1);

        if isempty(bIdx)
            continue
        end

        nTrials = C.all_nTrials_real(s, bIdx);
        xc = C.all_xc_real{s, bIdx};

        if isempty(nTrials) || isnan(nTrials) || nTrials < minTrials || isempty(xc)
            continue
        end

        lagsSec = C.lags(:) * C.binSize;
        allLagVals = [allLagVals; lagsSec(:)]; %#ok<AGROW>
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
        bIdx = find(C.behaviors == b, 1);

        if isempty(bIdx)
            continue
        end

        nTrials = C.all_nTrials_real(s, bIdx);
        if isempty(nTrials) || isnan(nTrials) || nTrials < minTrials
            continue
        end

        permLags = C.all_peakLagSec_perm{s, bIdx};
        if ~isempty(permLags)
            permLags = permLags(:);
            permLags = permLags(~isnan(permLags));
            allPermVals = [allPermVals; permLags]; %#ok<AGROW>
        end

        peakLag = C.all_peakLagSec_real(s, bIdx);
        if ~isempty(peakLag) && ~isnan(peakLag)
            allActualVals = [allActualVals; peakLag]; %#ok<AGROW>
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

%% ---- plot 1: per animal 1x3 lag vs correlation ----
for s = 1:nSess
    figure('Name', sprintf('%s chunked lag vs corr selected behaviors', animalIDs{s}), 'Color', 'w');
    tile_lay1 = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tile_lay1, sprintf('%s M1 Lag vs. Correlation', animalIDs{s}), ...
        'FontSize', titleFont, 'FontWeight', 'bold');

    firstLegendDone = false;

    for k = 1:numel(behNums)
        b = behNums(k);
        bIdx = find(C.behaviors == b, 1);

        nexttile(tile_lay1, k); hold on;

        if isempty(bIdx)
            title(behNamesPretty{k}, 'FontSize', panelTitleFont);
            axis off;
            continue
        end

        nTrials = C.all_nTrials_real(s, bIdx);
        xc = C.all_xc_real{s, bIdx};

        if isempty(nTrials) || isnan(nTrials) || nTrials < minTrials || isempty(xc)
            title(behNamesPretty{k}, 'FontSize', panelTitleFont);
            axis off;
            continue
        end

        permLags = C.all_peakLagSec_perm{s, bIdx};
        shiftZero = C.all_xcZeroLag_shift{s, bIdx};

        thisCorrCI = [NaN NaN];
        if ~isempty(shiftZero)
            shiftZero = shiftZero(~isnan(shiftZero));
            if ~isempty(shiftZero)
                thisCorrCI = prctile(shiftZero, [2.5 97.5]);
            end
        end

        thisLagCI = [NaN NaN];
        if ~isempty(permLags)
            permLags = permLags(~isnan(permLags));
            if ~isempty(permLags)
                thisLagCI = prctile(permLags, [2.5 97.5]);
            end
        end

        lagsSec = C.lags * C.binSize;
        peakLag = C.all_peakLagSec_real(s, bIdx);

        actualLagMat(s,k) = peakLag;
        if ~any(isnan(thisLagCI))
            lagCIMatLo(s,k) = thisLagCI(1);
            lagCIMatHi(s,k) = thisLagCI(2);
        end
        trialMat(s,k) = nTrials;

        hOrig = plot(lagsSec, xc, 'Color', origColor, 'LineWidth', 2);

        hCorr25 = gobjects(1);
        hCorr97 = gobjects(1);
        if ~any(isnan(thisCorrCI))
            hCorr25 = yline(thisCorrCI(1), '--', 'Color', corrCIColor, 'LineWidth', 1.4);
            hCorr97 = yline(thisCorrCI(2), '--', 'Color', corrCIColor, 'LineWidth', 1.4);
        end

        hPeak = xline(peakLag, '-', 'Color', peakLagColor, 'LineWidth', 1.8);

        hLag1 = gobjects(1);
        if ~any(isnan(thisLagCI))
            hLag1 = xline(thisLagCI(1), '--', 'Color', lagCIColor, 'LineWidth', 1.4);
            xline(thisLagCI(2), '--', 'Color', lagCIColor, 'LineWidth', 1.4);
        end

        xlim(lagXLimBySess(s,:));
        xlabel('Lag (Seconds)', 'FontSize', labelFont);
        ylabel('Correlation', 'FontSize', labelFont);
        title(behNamesPretty{k}, 'FontSize', panelTitleFont);
        box off;
        set(gca, 'FontSize', tickFont, 'LineWidth', 1, 'TickDir', 'out');

        if ~firstLegendDone
            legHandles = [hOrig];
            legLabels = {'Real Cross-Correlation'};
            if isgraphics(hCorr25), legHandles(end+1) = hCorr25; legLabels{end+1} = '2.5% Shift Control Correlation'; end
            if isgraphics(hCorr97), legHandles(end+1) = hCorr97; legLabels{end+1} = '97.5% Shift Control Correlation'; end
            legHandles(end+1) = hPeak; legLabels{end+1} = 'Actual Peak Lag';
            if isgraphics(hLag1), legHandles(end+1) = hLag1; legLabels{end+1} = '95% Permutation Lag Bounds'; end

            lgd = legend(legHandles, legLabels, 'Orientation', 'horizontal');
            lgd.Layout.Tile = 'south';
            lgd.FontSize = legendFont;
            lgd.Box = 'off';
            firstLegendDone = true;
        end
    end
end

%% ---- plot 3: per animal 1x3 permutation histograms ----
for s = 1:nSess
    figure('Name', sprintf('%s chunked permutation selected behaviors', animalIDs{s}), 'Color', 'w');
    tile_lay2 = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tile_lay2, sprintf('%s Peak Lags vs. Permutation Distribution', animalIDs{s}), ...
        'FontSize', titleFont, 'FontWeight', 'bold');

    sharedLegendHandles = gobjects(4,1);
    legendSet = false;

    for k = 1:numel(behNums)
        b = behNums(k);
        bIdx = find(C.behaviors == b, 1);

        nexttile(tile_lay2, k); hold on;

        if isempty(bIdx)
            title(behNamesPretty{k}, 'FontSize', panelTitleFont);
            axis off;
            continue
        end

        nTrials = C.all_nTrials_real(s, bIdx);
        permLags = C.all_peakLagSec_perm{s, bIdx};

        if isempty(nTrials) || isnan(nTrials) || nTrials < minTrials || isempty(permLags)
            title(behNamesPretty{k}, 'FontSize', panelTitleFont);
            axis off;
            continue
        end

        permLags = permLags(~isnan(permLags));
        if isempty(permLags)
            title(behNamesPretty{k}, 'FontSize', panelTitleFont);
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
        hPeak = xline(C.all_peakLagSec_real(s, bIdx), 'r-', 'LineWidth', 1.8);

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
        title(behNamesPretty{k}, 'FontSize', panelTitleFont);
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

%% ---- plot 2: one summary 1x3 figure with all animals in each behavior panel ----
figure('Name', 'Chunked selected behavior summary', 'Color', 'w');
tile_lay3 = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tile_lay3, 'Actual Peak Lags vs. Permutation Control Range', ...
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
    title(behNamesPretty{k}, 'FontSize', panelTitleFont);
    box off;
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
