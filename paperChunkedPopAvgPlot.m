function paperChunkedPopAvgPlot(saveFile)
% plots saved trial-averaged chunked population-averaged cross-correlation outputs (but no d043!)
% run: paperChunkedPopAvgPlot("C:\Users\mirilab\Documents\GlobusTransfer\concatCrossCorr_trialavg_chunked_popavg_ALL_SESSIONS.mat")

arguments
    saveFile (1,1) string
end

S = load(saveFile, 'allSessions');
allSessions = S.allSessions;
sessions = allSessions.sessions;
nSess = numel(sessions);

%% ---- export directory ----
saveDir = "X:\David\AnalysesData\InterneuronAnalyses\AA Paper Plots";

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

origColor = [0 0 0];
corrCIColor = [0 0.2 0.6];
peakLagColor = [0.95 0.45 0.35];
lagCIColor = [0 0 0];
permHistColor = [0.3 0.6 0.8];

% Plot formatting
legendFont = 10;
titleFont = 10;
labelFont = 10;
tickFont = 10;
axesLineWidth = 0.5;

animalIDs = strings(1, nSess);

for iDir = 1:nSess
    sess = sessions(iDir);

    if isfield(sess, 'mouseID') && ~isempty(sess.mouseID)
        animalIDs(iDir) = string(sess.mouseID);
    elseif isfield(sess, 'sessionTag') && ~isempty(sess.sessionTag)
        animalIDs(iDir) = string(sess.sessionTag);
    else
        animalID = regexp(sess.baseDir, 'D\d+', 'match', 'once');

        if isempty(animalID)
            animalID = sprintf('Animal %d', iDir);
        end

        animalIDs(iDir) = string(animalID);
    end
end

%% ---- remove D043 from plotting ----
keepSess = ~strcmpi(animalIDs, "D043");
sessions = sessions(keepSess);
animalIDs = animalIDs(keepSess);
nSess = numel(sessions);

peakLags = nan(nSess,1);
lagCIAll = nan(nSess,2);
permLagCell = cell(nSess,1);

%% ---- tiled cross-correlation figure ----
fig1 = figure( ...
    'Name', 'M1 lag vs. correlation', ...
    'Color', 'w');

tile_lay = tiledlayout(1, nSess, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

title(tile_lay, ...
    'M1 lag vs. correlation', ...
    'FontSize', titleFont);

for iDir = 1:nSess
    nexttile(tile_lay, iDir);
    hold on;

    sess = sessions(iDir);

    if isempty(sess.real_xc) || isempty(sess.lags)
        title(sprintf('%s (missing)', animalIDs(iDir)), ...
            'FontSize', titleFont);
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
            corrCI = prctile( ...
                sess.shift_xcZeroLag(goodShift), ...
                [2.5 97.5]);
        else
            corrCI = [NaN NaN];
        end
    else
        corrCI = [NaN NaN];
    end

    if ~isempty(sess.perm_peakLagSec)
        goodPerm = ~isnan(sess.perm_peakLagSec);

        if any(goodPerm)
            lagCI = prctile( ...
                sess.perm_peakLagSec(goodPerm), ...
                [2.5 97.5]);

            permLagCell{iDir} = ...
                sess.perm_peakLagSec(goodPerm);
        else
            lagCI = [NaN NaN];
        end
    else
        lagCI = [NaN NaN];
    end

    lagCIAll(iDir,:) = lagCI;

    hOrig = plot(lagsSec, xc, ...
        'Color', origColor, ...
        'LineWidth', 2);

    hCorr25 = gobjects(1);
    hCorr97 = gobjects(1);

    if ~any(isnan(corrCI))
        hCorr25 = yline(corrCI(1), '--', ...
            'Color', corrCIColor, ...
            'LineWidth', 1.4);

        hCorr97 = yline(corrCI(2), '--', ...
            'Color', corrCIColor, ...
            'LineWidth', 1.4);
    end

    hPeak = xline(peakLag, '-', ...
        'Color', peakLagColor, ...
        'LineWidth', 1.8);

    hLagCI = gobjects(1,2);

    if ~any(isnan(lagCI))
        hLagCI(1) = xline(lagCI(1), '--', ...
            'Color', lagCIColor, ...
            'LineWidth', 1.4);

        hLagCI(2) = xline(lagCI(2), '--', ...
            'Color', lagCIColor, ...
            'LineWidth', 1.4);
    end

    xlabel('Lag (seconds)', ...
        'FontSize', labelFont);

    ylabel('Correlation', ...
        'FontSize', labelFont);

    title(animalIDs(iDir), ...
        'FontSize', titleFont);

    box off;

    set(gca, ...
        'FontSize', tickFont, ...
        'LineWidth', axesLineWidth, ...
        'TickDir', 'out', ...
        'XColor', 'k', ...
        'YColor', 'k');

    if iDir == 1
        legendHandles = hOrig;
        legendLabels = {'Real cross-correlation'};

        if isgraphics(hCorr25)
            legendHandles(end+1) = hCorr25;
            legendLabels{end+1} = ...
                '2.5% shift control correlation';
        end

        if isgraphics(hCorr97)
            legendHandles(end+1) = hCorr97;
            legendLabels{end+1} = ...
                '97.5% shift control correlation';
        end

        legendHandles(end+1) = hPeak;
        legendLabels{end+1} = ...
            'Actual peak lag';

        if isgraphics(hLagCI(1))
            legendHandles(end+1) = hLagCI(1);
            legendLabels{end+1} = ...
                '95% permutation lag bounds';
        end

        lgd = legend( ...
            legendHandles, ...
            legendLabels, ...
            'Orientation', 'horizontal');

        lgd.Layout.Tile = 'south';
        lgd.FontSize = legendFont;
        lgd.Box = 'off';
    end
end

% Export entire tiled figure as vector PDF
exportgraphics(fig1, ...
    fullfile(saveDir, 'paper_cross_correlation.pdf'), ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');


%% ---- summary figure ----
fig2 = figure( ...
    'Name', 'Actual peak lags vs. permutation control range', ...
    'Color', 'w');

hold on;

xPos = 1:nSess;

hRange = gobjects(1);

for i = 1:nSess
    if ~any(isnan(lagCIAll(i,:)))
        hRange = line( ...
            [xPos(i) xPos(i)], ...
            lagCIAll(i,:), ...
            'Color', [0.6 0.6 0.6], ...
            'LineWidth', 2);
    end
end

hActual = scatter(xPos, peakLags, 70, 'k', 'filled');

yline(0, 'k:');

xlim([0.5 nSess + 0.5]);

xlabel('Animal', ...
    'FontSize', labelFont);

ylabel('Peak lag (s)', ...
    'FontSize', labelFont);

xticks(xPos);
xticklabels(cellstr(animalIDs));

title( ...
    'Actual peak lags vs. permutation control range', ...
    'FontSize', titleFont);

box off;
grid on;

set(gca, ...
    'FontSize', tickFont, ...
    'LineWidth', axesLineWidth, ...
    'TickDir', 'out', ...
    'XColor', 'k', ...
    'YColor', 'k');

if isgraphics(hRange)
    legend([hActual hRange], ...
        {'Actual peak lag', ...
         'Permutation control range'}, ...
        'Location', 'southoutside', ...
        'Orientation', 'horizontal', ...
        'Box', 'off', ...
        'FontSize', legendFont);
end

% Export summary figure as vector PDF
exportgraphics(fig2, ...
    fullfile(saveDir, 'paper_peak_lag_summary.pdf'), ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');


%% ---- common histogram limits ----
allPermLagsForLimits = cat(2, permLagCell{:});
allPermLagsForLimits = ...
    allPermLagsForLimits(~isnan(allPermLagsForLimits));

allActualPeakLags = peakLags(~isnan(peakLags));

if ~isempty(allPermLagsForLimits) || ...
        ~isempty(allActualPeakLags)

    combinedLagVals = ...
        [allPermLagsForLimits(:); allActualPeakLags(:)];

    xMinCommon = min(combinedLagVals);
    xMaxCommon = max(combinedLagVals);

    if xMinCommon == xMaxCommon
        xPad = 0.001;
    else
        xPad = 0.05 * ...
            (xMaxCommon - xMinCommon);
    end

    commonXLim = ...
        [xMinCommon - xPad, xMaxCommon + xPad];

else
    commonXLim = [-0.1 0.1];
end

nBins = 24;

commonEdges = linspace( ...
    commonXLim(1), ...
    commonXLim(2), ...
    nBins + 1);


%% ---- per-animal permutation histograms ----
fig3 = figure( ...
    'Name', 'Peak lags vs. permutation distribution', ...
    'Color', 'w');

tile_lay2 = tiledlayout( ...
    1, nSess, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

title(tile_lay2, ...
    'Peak lags vs. permutation distribution', ...
    'FontSize', titleFont);

sharedLegendHandles = gobjects(4,1);
legendSet = false;

for s = 1:nSess

    ax = nexttile;
    hold(ax, 'on');

    permLags = permLagCell{s};
    permLags = permLags(~isnan(permLags));

    if isempty(permLags)
        title( ...
            sprintf('%s (no perms)', animalIDs(s)), ...
            'FontSize', titleFont);

        axis off;
        continue;
    end

    hHist = histogram(permLags, ...
        'BinEdges', commonEdges, ...
        'FaceColor', permHistColor, ...
        'EdgeColor', 'none');

    prcLag = prctile( ...
        permLags, ...
        [2.5 97.5]);

    hPrc1 = xline(prcLag(1), '--', ...
        'Color', [0.2 0.2 0.2], ...
        'LineWidth', 1.5);

    hPrc2 = xline(prcLag(2), '--', ...
        'Color', [0.2 0.2 0.2], ...
        'LineWidth', 1.5);

    hActual = xline(peakLags(s), 'r-', ...
        'LineWidth', 1.8);

    if ~legendSet
        sharedLegendHandles(1) = hActual;

        sharedLegendHandles(2) = patch( ...
            nan, nan, permHistColor, ...
            'EdgeColor', 'none');

        sharedLegendHandles(3) = hPrc1;
        sharedLegendHandles(4) = hPrc2;

        legendSet = true;
    end

    xlim(commonXLim);

    xlabel('Peak lag (s)', ...
        'FontSize', labelFont);

    ylabel('Count', ...
        'FontSize', labelFont);

    title(animalIDs(s), ...
        'FontSize', titleFont);

    box off;

    set(gca, ...
        'FontSize', tickFont, ...
        'LineWidth', axesLineWidth, ...
        'TickDir', 'out', ...
        'XColor', 'k', ...
        'YColor', 'k');
end

if legendSet

    lgd2 = legend( ...
        sharedLegendHandles, ...
        {'Actual peak lag', ...
         'Permuted peak lags', ...
         '2.5% permutation control lag', ...
         '97.5% permutation control lag'}, ...
        'Orientation', 'horizontal');

    lgd2.Layout.Tile = 'south';
    lgd2.FontSize = legendFont;
    lgd2.Box = 'off';
end

% Export entire tiled histogram figure as vector PDF
exportgraphics(fig3, ...
    fullfile(saveDir, ...
    'paper_per_animal_permutation_distributions.pdf'), ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');


%% ---- combined histogram ----
allPermLags = cat(2, permLagCell{:});
allPermLags = allPermLags(~isnan(allPermLags));

fig4 = figure( ...
    'Name', 'All animals: peak lags vs. permutation distribution', ...
    'Color', 'w');

hold on;

hHist = histogram(allPermLags, ...
    'BinEdges', commonEdges, ...
    'FaceColor', permHistColor, ...
    'EdgeColor', 'none');

xlabel('Peak lag (s)', ...
    'FontSize', labelFont);

ylabel('Count', ...
    'FontSize', labelFont);

prcAll = prctile( ...
    allPermLags, ...
    [2.5 97.5]);

h2_5 = xline(prcAll(1), '--', ...
    'Color', [0.2 0.2 0.2], ...
    'LineWidth', 1.5);

h97_5 = xline(prcAll(2), '--', ...
    'Color', [0.2 0.2 0.2], ...
    'LineWidth', 1.5);

co = lines(nSess);
actualLines = gobjects(0,1);

for s = 1:nSess

    if ~isnan(peakLags(s))
        actualLines(end+1,1) = ...
            xline(peakLags(s), '-', ...
            'Color', co(s,:), ...
            'LineWidth', 1.5);
    end
end

xlim(commonXLim);

title( ...
    'All animals: peak lags vs. permutation distribution', ...
    'FontSize', titleFont);

box off;

set(gca, ...
    'FontSize', tickFont, ...
    'LineWidth', axesLineWidth, ...
    'TickDir', 'out', ...
    'XColor', 'k', ...
    'YColor', 'k');

legEntries = {
    'Permuted peak lags', ...
    '2.5% permutation control lag', ...
    '97.5% permutation control lag'
};

for s = 1:nSess
    legEntries{end+1} = ...
        sprintf('%s actual peak lag', animalIDs(s));
end

legend( ...
    [hHist h2_5 h97_5 actualLines(:)'], ...
    legEntries, ...
    'Location', 'best', ...
    'Box', 'off', ...
    'FontSize', legendFont);

% Export combined histogram as vector PDF
exportgraphics(fig4, ...
    fullfile(saveDir, ...
    'paper_combined_permutation_distribution.pdf'), ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');

end
