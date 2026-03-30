function plotSavedChunkedPopavgCrossCorrelationResults(saveFile)
% plots saved chunked population-averaged cross-correlation outputs without rerunning the analysis

% updated to match the plotting style of plotSavedCrossCorrelationResults that is specific to my non-chunked population-averaged results

% run: plotSavedChunkedPopavgCrossCorrelationResults("C:\Users\mirilab\Documents\GlobusTransfer\concatCrossCorr_chunked_popavg_ALL_SESSIONS.mat")

arguments
    saveFile (1,1) string
end

S = load(saveFile, 'allSessions');
allSessions = S.allSessions;

sessions = allSessions.sessions;
nSess = numel(sessions);

% colors
origColor = [0 0 0];
corrCIColor = [0 0.2 0.6];
peakLagColor = [0.95 0.45 0.35];
lagCIColor = [0 0 0];
permHistColor = [0.3 0.6 0.8];

% font sizes
legendFont = 17;
titleFont = 19;
labelFont = 19;
tickFont = 17;

peakLags = nan(nSess,1);
lagCIAll = nan(nSess,2);
permLagCell = cell(nSess,1);

%% ---- get animal ids for all sessions ----
animalIDs = strings(1, nSess);
for iDir = 1:nSess
    sess = sessions(iDir);
    animalID = regexp(sess.baseDir, 'D\d+', 'match', 'once');
    if isempty(animalID)
        animalID = sprintf('Animal %d', iDir);
    end
    animalIDs(iDir) = string(animalID);
end

%% ---- tiled cross-correlation figure ----
figure('Name', 'M1 Lag vs. Correlation', 'Color', 'w');
tile_lay = tiledlayout(1, nSess, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tile_lay, 'M1 Lag vs. Correlation', 'FontSize', titleFont);

for iDir = 1:nSess
    nexttile(tile_lay, iDir);
    hold on;

    sess = sessions(iDir);

    if isempty(sess.real_xc) || isempty(sess.lags)
        title(sprintf('%s (missing)', animalIDs(iDir)), 'FontSize', titleFont);
        axis off;
        continue;
    end

    lagsSec = sess.lags * sess.binSize;
    xc = sess.real_xc(:)';
    peakLag = sess.real_peakLagSec;
    peakLags(iDir) = peakLag;

    if ~isempty(sess.shift_xcZeroLag)
        goodShift = ~isnan(sess.shift_xcZeroLag);
        if any(goodShift)
            corrCI = prctile(sess.shift_xcZeroLag(goodShift), [2.5 97.5]);
        else
            corrCI = [NaN NaN];
        end
    else
        corrCI = [NaN NaN];
    end

    if ~isempty(sess.perm_peakLagSec)
        goodPerm = ~isnan(sess.perm_peakLagSec);
        if any(goodPerm)
            lagCI = prctile(sess.perm_peakLagSec(goodPerm), [2.5 97.5]);
            permLagCell{iDir} = sess.perm_peakLagSec(goodPerm);
        else
            lagCI = [NaN NaN];
        end
    else
        lagCI = [NaN NaN];
    end
    lagCIAll(iDir,:) = lagCI;

    hOrig = plot(lagsSec, xc, 'Color', origColor, 'LineWidth', 2);

    hCorr25 = gobjects(1);
    hCorr97 = gobjects(1);
    if ~any(isnan(corrCI))
        hCorr25 = yline(corrCI(1), '--', 'Color', corrCIColor, 'LineWidth', 1.4);
        hCorr97 = yline(corrCI(2), '--', 'Color', corrCIColor, 'LineWidth', 1.4);
    end

    hPeak = xline(peakLag, '--', 'Color', peakLagColor, 'LineWidth', 1.8);

    hLagCI = gobjects(1,2);
    if ~any(isnan(lagCI))
        hLagCI(1) = xline(lagCI(1), '--', 'Color', lagCIColor, 'LineWidth', 1.4);
        hLagCI(2) = xline(lagCI(2), '--', 'Color', lagCIColor, 'LineWidth', 1.4);
    end

    xlabel('Lag (Seconds)', 'FontSize', labelFont);
    ylabel('Correlation', 'FontSize', labelFont);
    title(animalIDs(iDir), 'FontSize', titleFont);
    box off;
    set(gca, 'FontSize', tickFont, 'LineWidth', 1, 'TickDir', 'out');

    if iDir == 1
        legendHandles = hOrig;
        legendLabels = {'Real Cross-Correlation'};

        if isgraphics(hCorr25)
            legendHandles(end+1) = hCorr25; %#ok<AGROW>
            legendLabels{end+1} = '2.5% Shift Control Correlation'; %#ok<AGROW>
        end
        if isgraphics(hCorr97)
            legendHandles(end+1) = hCorr97; %#ok<AGROW>
            legendLabels{end+1} = '97.5% Shift Control Correlation'; %#ok<AGROW>
        end

        legendHandles(end+1) = hPeak; %#ok<AGROW>
        legendLabels{end+1} = 'Actual Peak Lag'; %#ok<AGROW>

        if isgraphics(hLagCI(1))
            legendHandles(end+1) = hLagCI(1); %#ok<AGROW>
            legendLabels{end+1} = 'Lag 95% CI (Permutation)'; %#ok<AGROW>
        end

        lgd = legend(legendHandles, legendLabels, 'Orientation', 'horizontal');
        lgd.Layout.Tile = 'south';
        lgd.FontSize = legendFont;
        lgd.Box = 'off';
    end
end

%% ---- summary figure: peak lag + ci per animal ----
figure('Name','Cortex peak lag summary with permutation CI', 'Color', 'w');
hold on;

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

nBins = 24;
commonEdges = linspace(commonXLim(1), commonXLim(2), nBins + 1);

%% ---- per-animal histograms ----
figure('Name','Cortex chunked XC peak lag permutation distributions (per animal)', ...
       'Color','w');
tile_lay2 = tiledlayout(1, nSess, 'TileSpacing','compact','Padding','compact');
title(tile_lay2, 'Actual Peak Lags vs. Permutation Distribution', 'FontSize', titleFont);

sharedLegendHandles = gobjects(4,1);
legendSet = false;

for s = 1:nSess
    ax = nexttile;
    hold(ax, 'on');

    permLags = permLagCell{s};
    permLags = permLags(~isnan(permLags));

    if isempty(permLags)
        title(sprintf('%s (no perms)', animalIDs(s)), 'FontSize', titleFont);
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

    if ~legendSet
        sharedLegendHandles(1) = hActual;

        hHistLegend = patch(nan, nan, permHistColor, 'EdgeColor', 'none');
        sharedLegendHandles(2) = hHistLegend;

        sharedLegendHandles(3) = hPrc1;
        sharedLegendHandles(4) = hPrc2;

        legendSet = true;
    end

    xlim(commonXLim);
    xlabel('Peak Lag (s)', 'FontSize', labelFont);
    ylabel('Count', 'FontSize', labelFont);
    title(animalIDs(s), 'FontSize', titleFont);
    box off;
    set(gca, 'FontSize', tickFont, 'LineWidth', 1, 'TickDir', 'out');
end

if legendSet
    lgd2 = legend(sharedLegendHandles, ...
        {'Actual Peak Lag', ...
         'Permuted Peak Lags', ...
         '2.5% Permutation Control Lag', ...
         '97.5% Permutation Control Lag'}, ...
        'Orientation', 'horizontal');
    lgd2.Layout.Tile = 'south';
    lgd2.FontSize = legendFont;
    lgd2.Box = 'off';
end

%% ---- combined histogram across animals ----
allPermLags = cat(2, permLagCell{:});
validAll = ~isnan(allPermLags);
allPermLags = allPermLags(validAll);

figure('Name','Cortex chunked XC peak lag permutations (all animals combined)', ...
       'Color','w');
hold on;

hHist = histogram(allPermLags, ...
    'BinEdges', commonEdges, ...
    'FaceColor', permHistColor, ...
    'EdgeColor','none');

xlabel('Peak Lag (s)', 'FontSize', labelFont);
ylabel('Count', 'FontSize', labelFont);

prcAll = prctile(allPermLags, [2.5 97.5]);
h2_5 = xline(prcAll(1), '--', 'Color',[0.2 0.2 0.2], 'LineWidth',1.5);
h97_5 = xline(prcAll(2), '--', 'Color',[0.2 0.2 0.2], 'LineWidth',1.5);

co = lines(nSess);
actualLines = gobjects(0,1);

for s = 1:nSess
    if ~isnan(peakLags(s))
        actualLines(end+1,1) = xline(peakLags(s), '-', ...
            'Color',co(s,:), 'LineWidth',1.5); %#ok<AGROW>
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
