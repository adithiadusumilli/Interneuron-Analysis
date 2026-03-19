function plotSavedCrossCorrelationResults(saveFile)
% plotsavedcrosscorrelationresults for non-chunked xc

% loads saved pop-averaged cross-correlation outputs and redraws figures without rerunning analysis

% edits in this version:
%   - uses animal names in titles instead of session numbers
%   - adds one shared legend below the per-animal permutation histogram figure
%   - forces the same x-axis limits across all per-animal permutation histogram panels
%   - reduces histogram bin size
%   - increases font size across figures
%   - legend text for percentile lines does not include numeric values

% run:
% plotSavedCrossCorrelationResults("X:\David\AnalysesData\InterneuronAnalyses\Lab Meeting Pres\4 aninals run cross correlation, no chunking, pop-wise\runCrossCorrelation_savedOutputs_all4Animals.mat")

arguments
    saveFile (1,1) string
end

S = load(saveFile, 'xcorrResults', 'peakLags', 'peakCorrs', 'lagCIAll', 'permLagCell', 'animalLabels');

xcorrResults = S.xcorrResults;
peakLags = S.peakLags;
peakCorrs = S.peakCorrs; %#ok<NASGU>
lagCIAll = S.lagCIAll;
permLagCell = S.permLagCell;

nSess = numel(xcorrResults.sessions);

% colors
origColor = [0 0 0];
corrCIColor = [0 0.2 0.6];
peakLagColor = [0.95 0.45 0.35];
lagCIColor = [0 0 0];
permHistColor = [0.3 0.6 0.8];

baseFont = 20;
legendFont = 18;
titleFont = 20;
labelFont = 20;
tickFont = 18;

%% ---- get animal ids for all sessions ----
animalIDs = strings(1, nSess);
for iDir = 1:nSess
    sess = xcorrResults.sessions(iDir);
    animalID = regexp(sess.baseDir, 'D\d+', 'match', 'once');
    if isempty(animalID)
        animalID = sprintf('Animal %d', iDir);
    end
    animalIDs(iDir) = string(animalID);
end

%% ---- tiled cross-correlation figure ----
figure('Name', 'M1 Lag vs. Correlation', 'Color', 'w');
tile_lay = tiledlayout(1, nSess, 'TileSpacing', 'compact', 'Padding', 'compact');

for iDir = 1:nSess
    nexttile(tile_lay, iDir); hold on;

    sess = xcorrResults.sessions(iDir);

    if isempty(sess.lagsSec) || isempty(sess.xc)
        title(sprintf('%s (missing)', animalIDs(iDir)), 'FontSize', titleFont);
        axis off;
        continue;
    end

    animalID = animalIDs(iDir);

    hOrig = plot(sess.lagsSec, sess.xc, 'Color', origColor, 'LineWidth', 2);
    hCorr25 = yline(sess.corrCI(1), '--', 'Color', corrCIColor, 'LineWidth', 1.4);
    hCorr97 = yline(sess.corrCI(2), '--', 'Color', corrCIColor, 'LineWidth', 1.4);
    hPeak = xline(sess.peakLag, '--', 'Color', peakLagColor, 'LineWidth', 1.8);

    hLagCI = gobjects(1,2);
    if ~any(isnan(sess.lagCI))
        hLagCI(1) = xline(sess.lagCI(1), '--', 'Color', lagCIColor, 'LineWidth', 1.4);
        hLagCI(2) = xline(sess.lagCI(2), '--', 'Color', lagCIColor, 'LineWidth', 1.4);
    end

    xlabel('Lag (Seconds)', 'FontSize', labelFont);
    ylabel('Correlation', 'FontSize', labelFont);
    title(sprintf('M1 Lag vs. Correlation — %s', animalID), 'FontSize', titleFont);
    box off;
    set(gca, 'FontSize', tickFont, 'LineWidth', 1, 'TickDir', 'out');

    if iDir == 1
        lgd = legend([hOrig hCorr25 hCorr97 hPeak hLagCI(1)], ...
            {'Real Cross-Correlation', ...
             '2.5% Shift Control Correlation', ...
             '97.5% Shift Control Correlation', ...
             'Actual Peak Lag', ...
             'Lag 95% CI (Permutation)'}, ...
            'Orientation', 'horizontal');
        lgd.Layout.Tile = 'south';
        lgd.FontSize = legendFont;
        lgd.Box = 'off';
    end
end

%% ---- summary figure: peak lag + ci per animal ----
figure('Name','Cortex peak lag summary with permutation CI', 'Color', 'w'); hold on;

xPos = 1:nSess;

for i = 1:nSess
    if ~any(isnan(lagCIAll(i,:)))
        line([xPos(i) xPos(i)], lagCIAll(i,:), 'Color', [0.6 0.6 0.6], 'LineWidth', 2);
    end
end

scatter(xPos, peakLags, 70, 'k', 'filled');
yline(0, 'k:');

xlim([0.5 nSess + 0.5]);
xlabel('Animal', 'FontSize', labelFont);
ylabel('Peak Lag (Seconds)', 'FontSize', labelFont);
xticks(xPos);
xticklabels(cellstr(animalIDs));

title('Cortex Peak Lags with 95% Permutation CI', 'FontSize', titleFont);
box off;
grid on;
set(gca, 'FontSize', tickFont, 'LineWidth', 1, 'TickDir', 'out');

%% ---- set common histogram x-limits across per-animal panels ----
allPermLagsForLimits = cat(2, permLagCell{:});
allPermLagsForLimits = allPermLagsForLimits(~isnan(allPermLagsForLimits));

allActualPeakLags = peakLags(~isnan(peakLags));

if ~isempty(allPermLagsForLimits) || ~isempty(allActualPeakLags)
    combinedLagVals = [allPermLagsForLimits(:); allActualPeakLags(:)];
    xMinCommon = min(combinedLagVals);
    xMaxCommon = max(combinedLagVals);

    if xMinCommon == xMaxCommon
        xPad = 0.001;
    else
        xPad = 0.05 * (xMaxCommon - xMinCommon);
    end

    commonXLim = [xMinCommon - xPad, xMaxCommon + xPad];
else
    commonXLim = [-0.1 0.1];
end

% reduce bin size
nBins = 24;
commonEdges = linspace(commonXLim(1), commonXLim(2), nBins + 1);

%% ---- per-animal histograms ----
figure('Name','Cortex non-chunked XC peak lag permutation distributions (per animal)', ...
       'Color','w');
tile_lay2 = tiledlayout(1, nSess, 'TileSpacing','compact','Padding','compact');

sharedLegendHandles = gobjects(4,1);

for s = 1:nSess
    ax = nexttile; hold(ax, 'on');

    permLags = permLagCell{s};
    permLags = permLags(~isnan(permLags));

    if isempty(permLags)
        title(sprintf('%s: Actual Peak Lag vs. Permutations (no perms)', animalIDs(s)), ...
            'FontSize', titleFont);
        axis off;
        continue;
    end

    hHist = histogram(permLags, ...
        'BinEdges', commonEdges, ...
        'FaceColor', permHistColor, ...
        'EdgeColor', 'none');

    prcLag = prctile(permLags, [2.5 97.5]);
    hPrc1 = xline(prcLag(1), '--', 'Color',[0.2 0.2 0.2], 'LineWidth',1.5);
    hPrc2 = xline(prcLag(2), '--', 'Color',[0.2 0.2 0.2], 'LineWidth',1.5);
    hActual = xline(peakLags(s), 'r-', 'LineWidth',1.8);

    % store shared legend handles from first valid panel only
    if s == 1
        sharedLegendHandles(1) = hActual;
        sharedLegendHandles(2) = hHist;
        sharedLegendHandles(3) = hPrc1;
        sharedLegendHandles(4) = hPrc2;
    end

    xlim(commonXLim);
    xlabel('Peak Lag (s)', 'FontSize', labelFont);
    ylabel('Count', 'FontSize', labelFont);
    title(sprintf('%s: Actual Peak Lag vs. Permutations', animalIDs(s)), ...
        'FontSize', titleFont);
    box off;
    set(gca, 'FontSize', tickFont, 'LineWidth', 1, 'TickDir', 'out');
end

lgd2 = legend(sharedLegendHandles, ...
    {'Actual Peak Lag', ...
     'Permuted Peak Lags', ...
     '2.5% Permutation Control Lag', ...
     '97.5% Permutation Control Lag'}, ...
    'Orientation', 'horizontal');
lgd2.Layout.Tile = 'south';
lgd2.FontSize = legendFont;
lgd2.Box = 'off';

%% ---- combined histogram across animals ----
allPermLags = cat(2, permLagCell{:});
validAll    = ~isnan(allPermLags);
allPermLags = allPermLags(validAll);

figure('Name','Cortex non-chunked XC peak lag permutations (all animals combined)', ...
       'Color','w'); hold on;

hHist = histogram(allPermLags, ...
    'BinEdges', commonEdges, ...
    'FaceColor', permHistColor, ...
    'EdgeColor','none');
xlabel('Peak Lag (s)', 'FontSize', labelFont);
ylabel('Count', 'FontSize', labelFont);

prcAll = prctile(allPermLags, [2.5 97.5]);
h2_5   = xline(prcAll(1), '--', 'Color',[0.2 0.2 0.2], 'LineWidth',1.5);
h97_5  = xline(prcAll(2), '--', 'Color',[0.2 0.2 0.2], 'LineWidth',1.5);

co          = lines(nSess);
actualLines = gobjects(0,1);
for s = 1:nSess
    if ~isnan(peakLags(s))
        actualLines(end+1,1) = xline(peakLags(s), '-', 'Color',co(s,:), 'LineWidth',1.5); %#ok<AGROW>
    end
end

xlim(commonXLim);
title('All Animals Combined: Actual Peak Lags vs. Permutation Distribution', ...
    'FontSize', titleFont);
box off;
set(gca, 'FontSize', tickFont, 'LineWidth', 1, 'TickDir', 'out');

legEntries = {
    'Permuted Peak Lags', ...
    '2.5% Permutation Control Lag', ...
    '97.5% Permutation Control Lag'
};

for s = 1:nSess
    legEntries{end+1} = sprintf('%s Actual Peak Lag', animalIDs(s));
end

legend([hHist h2_5 h97_5 actualLines(:)'], legEntries, ...
    'Location', 'best', 'Box', 'off', 'FontSize', legendFont);

end
