function paperNoChunkPopAvgPlot(saveFile)
% plots saved non-chunked population-averaged cross-correlation outputs
% excludes D043 and exports paper-formatted vector PDFs

% run:
% paperNoChunkPopAvgPlot("X:\David\AnalysesData\InterneuronAnalyses\Lab Meeting Pres\6 animals run cross correlation, no chunking, pop-wise\runCrossCorrelation_savedOutputs_all6Animals.mat")

if nargin < 1 || isempty(saveFile)
    saveFile = ...
        "X:\David\AnalysesData\InterneuronAnalyses\Lab Meeting Pres\6 animals run cross correlation, no chunking, pop-wise\runCrossCorrelation_savedOutputs_all6Animals.mat";
end

S = load(saveFile);

xcorrResults = S.xcorrResults;
peakLags = S.peakLags;
peakCorrs = S.peakCorrs;
lagCIAll = S.lagCIAll;
permLagCell = S.permLagCell;

nSess = numel(peakLags);

if isfield(S, 'animalLabels')
    animalLabels = S.animalLabels;
else
    animalLabels = cell(1,nSess);
    for s = 1:nSess
        animalLabels{s} = sprintf('Session %d', s);
    end
end

if iscell(animalLabels)
    animalLabels = string(animalLabels);
end

%% ============================
%  remove D043
% ============================

keepSess = ~strcmpi(animalLabels, "D043");

animalLabels = animalLabels(keepSess);
peakLags = peakLags(keepSess);
peakCorrs = peakCorrs(keepSess);
lagCIAll = lagCIAll(keepSess,:);
permLagCell = permLagCell(keepSess);

% xcorrResults.sessions must be filtered as well
xcorrResults.sessions = xcorrResults.sessions(keepSess);

nSess = numel(peakLags);


%% ============================
%  export directory
% ============================

saveDir = ...
    "X:\David\AnalysesData\InterneuronAnalyses\AA Paper Plots\no-chunk popavg";

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end


%% ============================
%  plot formatting
% ============================

origColor = [0 0 0];
corrCIColor = [0.2 0.2 0.2];
peakLagColor = [0.95 0.45 0.35];
permHistColor = [0.3 0.6 0.8];
lagCIColor = [0.6 0.6 0.6];

legendFont = 10;
titleFont = 10;
labelFont = 10;
tickFont = 10;
axesLineWidth = 0.5;


%% ============================
%  plot cross-correlation curves
% ============================

fig1 = figure( ...
    'Name', 'M1 lag vs. correlation', ...
    'Color', 'w');

tile_lay = tiledlayout(1, nSess, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

title(tile_lay, ...
    'M1 lag vs. correlation', ...
    'FontSize', titleFont);

legendSet = false;

for s = 1:nSess

    nexttile(tile_lay, s);
    hold on;

    sess = xcorrResults.sessions(s);

    if isempty(sess.xc) || isempty(sess.lagsSec)
        title(sprintf('%s (missing)', animalLabels(s)), ...
            'FontSize', titleFont);
        axis off;
        continue;
    end

    hOrig = plot(sess.lagsSec, sess.xc, ...
        'Color', origColor, ...
        'LineWidth', 2);

    hPeak = xline(sess.peakLag, '-', ...
        'Color', peakLagColor, ...
        'LineWidth', 1.5);

    hCorr25 = yline(sess.corrCI(1), '--', ...
        'Color', corrCIColor, ...
        'LineWidth', 1);

    hCorr97 = yline(sess.corrCI(2), '--', ...
        'Color', corrCIColor, ...
        'LineWidth', 1);

    xlabel('Lag (seconds)', ...
        'FontSize', labelFont);

    ylabel('Correlation', ...
        'FontSize', labelFont);

    title(animalLabels(s), ...
        'FontSize', titleFont);

    box off;

    set(gca, ...
        'FontSize', tickFont, ...
        'LineWidth', axesLineWidth, ...
        'TickDir', 'out', ...
        'XColor', 'k', ...
        'YColor', 'k');

    if ~legendSet

        legendHandles = [ ...
            hOrig, ...
            hCorr25, ...
            hCorr97, ...
            hPeak];

        legendLabels = { ...
            'Real cross-correlation', ...
            '2.5% shift control correlation', ...
            '97.5% shift control correlation', ...
            'Actual peak lag'};

        lgd = legend( ...
            legendHandles, ...
            legendLabels, ...
            'Orientation', 'horizontal');

        lgd.Layout.Tile = 'south';
        lgd.FontSize = legendFont;
        lgd.Box = 'off';

        legendSet = true;
    end
end

% Export as vector PDF
exportgraphics(fig1, ...
    fullfile(saveDir, 'paper_cross_correlation.pdf'), ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');


%% ============================
%  summary peak lag + CI
% ============================

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
            'Color', lagCIColor, ...
            'LineWidth', 2);
    end
end

hActual = scatter( ...
    xPos, ...
    peakLags, ...
    70, ...
    'k', ...
    'filled');

% Zero-lag reference line
yline(0, 'k:');

xlim([0.5 nSess + 0.5]);

xlabel('Animal', ...
    'FontSize', labelFont);

ylabel('Peak lag (s)', ...
    'FontSize', labelFont);

xticks(xPos);
xticklabels(cellstr(animalLabels));

title( ...
    'Actual peak lags vs. permutation control range', ...
    'FontSize', titleFont);

box off;

% No grid lines
grid off;

set(gca, ...
    'FontSize', tickFont, ...
    'LineWidth', axesLineWidth, ...
    'TickDir', 'out', ...
    'XColor', 'k', ...
    'YColor', 'k');

if isgraphics(hRange)

    legend( ...
        [hActual hRange], ...
        {'Actual peak lag', ...
         'Permutation control range'}, ...
        'Location', 'southoutside', ...
        'Orientation', 'horizontal', ...
        'Box', 'off', ...
        'FontSize', legendFont);
end

fprintf('\n========== summary cortex ==========\n');

fprintf('peak lags (s):\n');
disp(peakLags);

fprintf('peak corrs:\n');
disp(peakCorrs);

fprintf('lag CIs (s):\n');
disp(lagCIAll);

% Export as vector PDF
exportgraphics(fig2, ...
    fullfile(saveDir, 'paper_peak_lag_summary.pdf'), ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');


%% ============================
%  common histogram limits
% ============================

allPermForLimits = cat(2, permLagCell{:});
allPermForLimits = ...
    allPermForLimits(~isnan(allPermForLimits));

allActualForLimits = ...
    peakLags(~isnan(peakLags));

if isempty(allPermForLimits) && ...
        isempty(allActualForLimits)

    commonXLim = [-0.5 0.5];

else

    allVals = [ ...
        allPermForLimits(:); ...
        allActualForLimits(:)];

    xMin = min(allVals);
    xMax = max(allVals);

    if xMin == xMax
        pad = 0.001;
    else
        pad = 0.05 * (xMax - xMin);
    end

    commonXLim = ...
        [xMin - pad, xMax + pad];
end

nBins = 24;

commonEdges = linspace( ...
    commonXLim(1), ...
    commonXLim(2), ...
    nBins + 1);


%% ============================
%  per-animal permutation histograms
% ============================

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
            sprintf('%s (no perms)', animalLabels(s)), ...
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
        'Color', corrCIColor, ...
        'LineWidth', 1.5);

    hPrc2 = xline(prcLag(2), '--', ...
        'Color', corrCIColor, ...
        'LineWidth', 1.5);

    hActual = xline(peakLags(s), '-', ...
        'Color', peakLagColor, ...
        'LineWidth', 1.5);

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

    title(animalLabels(s), ...
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

% Export as vector PDF
exportgraphics(fig3, ...
    fullfile(saveDir, ...
    'paper_per_animal_permutation_distributions.pdf'), ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');


%% ============================
%  combined histogram
% ============================

allPermLags = cat(2, permLagCell{:});
allPermLags = ...
    allPermLags(~isnan(allPermLags));

fig4 = figure( ...
    'Name', ...
    'All animals: peak lags vs. permutation distribution', ...
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
    'Color', corrCIColor, ...
    'LineWidth', 1.5);

h97_5 = xline(prcAll(2), '--', ...
    'Color', corrCIColor, ...
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

legHandles = ...
    [hHist h2_5 h97_5 actualLines(:)'];

legEntries = { ...
    'Permuted peak lags', ...
    '2.5% permutation control lag', ...
    '97.5% permutation control lag'};

for s = 1:nSess

    if ~isnan(peakLags(s))
        legEntries{end+1} = ...
            sprintf('%s actual peak lag', ...
            animalLabels(s));
    end
end

legend( ...
    legHandles, ...
    legEntries, ...
    'Location', 'best', ...
    'Box', 'off', ...
    'FontSize', legendFont);

% Export as vector PDF
exportgraphics(fig4, ...
    fullfile(saveDir, ...
    'paper_combined_permutation_distribution.pdf'), ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');

end
