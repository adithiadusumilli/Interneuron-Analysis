function plotSavedIntPyrLagSkewMAFDR_NoChunk(resultsFile)
% plots saved no-chunk pairwise significance results from runIntPyrLagSkewMAFDRPermutationTest

% expected input:
%   resultsFile = output .mat from runIntPyrLagSkewMAFDRPermutationTest
%   for analysisType = "nochunk"

% figures per session:
%   1. significant pairwise peak lag map
%   2. distribution of significant peak lags
%   3. peak lag versus peak correlation
%   4. skew relative to null distribution

% summary figure:
%   - actual skew versus null 95% interval across sessions

% j run: plotSavedIntPyrLagSkewMAFDR_NoChunk("C:\Users\mirilab\Documents\GlobusTransfer\pairwise_nochunk_allPairs_ALL_SESSIONS_COMBINED_intPyrSkewMAFDR_nochunk.mat")

arguments
    resultsFile (1,1) string
end

S = load(resultsFile, 'results');
R = S.results;

if ~isfield(R, 'sessions') || isempty(R.sessions)
    error('results file does not contain a valid results.sessions field.');
end

if isfield(R, 'analysisType')
    if ~strcmpi(string(R.analysisType), "nochunk")
        error('this plotting function is for no-chunk saved outputs only.');
    end
end

nSess = numel(R.sessions);

summaryActualSkew = nan(nSess,1);
summaryNullCI = nan(nSess,2);
summaryAnimalID = cell(1,nSess);
summaryNSig = nan(nSess,1);
summaryNTotal = nan(nSess,1);

for s = 1:nSess
    sess = R.sessions{s};

    animalID = sess.animalID;
    if isempty(animalID)
        animalID = sprintf('Session%d', s);
    end
    summaryAnimalID{s} = animalID;

    actual = sess.actual;

    if isempty(actual) || ~isfield(actual, 'sigFDRMask')
        fprintf('session %d: missing actual significance results, skipping\n', s);
        continue;
    end

    nInt = sess.nInt;
    nPyr = sess.nPyr;

    % ================= HEATMAP =================
    sigMat = nan(nInt, nPyr);

    for k = 1:numel(actual.rows)
        r = actual.rows(k);
        c = actual.cols(k) - nInt;

        if c >= 1 && c <= nPyr && actual.sigFDRMask(k)
            sigMat(r,c) = actual.lagVals(k);
        end
    end

    figure('Name', sprintf('%s Significant Pairwise Peak Lag Map', animalID), 'Color', 'w');
    hImg = imagesc(sigMat);
    set(hImg, 'AlphaData', ~isnan(sigMat));
    set(gca, 'Color', 'k');
    axis xy;
    colormap(parula);

    cb = colorbar;
    ylabel(cb, 'Peak Lag (Seconds)', 'FontSize', 14);

    xlabel('Pyramidal Neuron Index', 'FontSize', 14);
    ylabel('Interneuron Index', 'FontSize', 14);

    title(sprintf('%s Significant Pairwise Peak Lag Map (%d / %d Pairs)', ...
        animalID, actual.nSigFDR, actual.nPairsNominal), ...
        'FontSize', 16, 'FontWeight', 'bold');

    set(gca, 'FontSize', 13, 'TickDir', 'out');
    box off;

    % ================= LAG HIST =================
    figure('Name', sprintf('%s Distribution Of Peak Lags', animalID), 'Color', 'w');
    if ~isempty(actual.sigLagVec)
        edges = makeLagEdges(actual.sigLagVec);
        histogram(actual.sigLagVec, 'BinEdges', edges, 'EdgeColor', 'none');
    else
        histogram(nan);
    end
    hold on;

    hZero = xline(0, 'k--', 'LineWidth', 1.5);

    xlabel('Peak Lag (Seconds)', 'FontSize', 14);
    ylabel('Pair Count', 'FontSize', 14);

    title(sprintf('%s Distribution Of Peak Lags (%d / %d Significant Pairs)', ...
        animalID, actual.nSigFDR, actual.nPairsNominal), ...
        'FontSize', 16, 'FontWeight', 'bold');

    legend(hZero, 'Zero Lag Reference', 'Location', 'best', 'FontSize', 12);

    set(gca, 'FontSize', 13);
    grid on;

    % ================= SCATTER =================
    if ~isempty(actual.sigLagVec) && ~isempty(actual.realVals)
        sigCorrVec = actual.realVals(actual.sigFDRMask);
        sigCorrVec = sigCorrVec(~isnan(sigCorrVec) & isfinite(sigCorrVec));

        if numel(sigCorrVec) == numel(actual.sigLagVec)
            figure('Name', sprintf('%s Peak Lag Versus Peak Correlation', animalID), 'Color', 'w');
            scatterhist(actual.sigLagVec(:), sigCorrVec(:), ...
                'Direction', 'out', 'Marker', '.');

            hold on;

            hThresh = gobjects(1);
            if isfield(R, 'corrThresh') && ~isempty(R.corrThresh) && ~isnan(R.corrThresh)
                hThresh = yline(R.corrThresh, 'r--', 'LineWidth', 1.5);
            end

            xlabel('Peak Lag (Seconds)', 'FontSize', 14);
            ylabel('Peak Correlation', 'FontSize', 14);

            title(sprintf('%s Peak Lag Versus Peak Correlation (%d / %d Significant Pairs)', ...
                animalID, actual.nSigFDR, actual.nPairsNominal), ...
                'FontSize', 16, 'FontWeight', 'bold');

            if isgraphics(hThresh)
                legend(hThresh, 'Correlation Threshold', 'FontSize', 12, 'Location', 'best');
            end

            set(gca, 'FontSize', 13);
        end
    end

    % ================= SKEW =================
    validNullSkews = sess.validNullSkews;

    if isempty(validNullSkews)
        validNullSkews = nan(0,1);
    end

    if ~isempty(validNullSkews)
        nullCI = prctile(validNullSkews, [2.5 97.5]);
    else
        nullCI = [NaN NaN];
    end

    summaryActualSkew(s) = actual.skew;
    summaryNullCI(s,:) = nullCI;
    summaryNSig(s) = actual.nSigFDR;
    summaryNTotal(s) = actual.nPairsNominal;

    figure('Name', sprintf('%s Skew Relative To Null Distribution', animalID), 'Color', 'w');

    if ~isempty(validNullSkews)
        histogram(validNullSkews, 30, 'EdgeColor', 'none');
        hold on;
        hActual = xline(actual.skew, 'r', 'LineWidth', 2);
        hCI = xline(nullCI(1), 'k--', 'LineWidth', 1.5);
        xline(nullCI(2), 'k--', 'LineWidth', 1.5);

        legend([hActual hCI], ...
            {'Actual Skew', '95% Null Skew Interval'}, ...
            'FontSize', 12, 'Location', 'best');
    else
        histogram(nan);
        hold on;
        xline(actual.skew, 'r', 'LineWidth', 2);
    end

    xlabel('Skew (Mean Minus Median Over Standard Deviation)', 'FontSize', 14);
    ylabel('Count', 'FontSize', 14);

    title(sprintf('%s Skew Relative To Null Distribution', animalID), ...
        'FontSize', 16, 'FontWeight', 'bold');

    set(gca, 'FontSize', 13);
    grid on;
end

% ================= SUMMARY =================
figure('Name', 'Summary Of Skew Relative To Null Distribution Across Sessions', 'Color', 'w');
t = tiledlayout(1, nSess, 'TileSpacing', 'compact', 'Padding', 'compact');

for s = 1:nSess
    nexttile;
    hold on;

    xNull = 0.95;
    xActual = 1.05;

    hCI = gobjects(1);
    if ~any(isnan(summaryNullCI(s,:)))
        hCI = line([xNull xNull], summaryNullCI(s,:), ...
            'Color', [0.5 0.5 0.5], 'LineWidth', 4);

        plot(xNull, summaryNullCI(s,1), '_', 'Color', [0.3 0.3 0.3], ...
            'MarkerSize', 18, 'LineWidth', 2);
        plot(xNull, summaryNullCI(s,2), '_', 'Color', [0.3 0.3 0.3], ...
            'MarkerSize', 18, 'LineWidth', 2);
    end

    hDot = plot(xActual, summaryActualSkew(s), ...
        'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 9);

    yline(0, 'k:', 'LineWidth', 1.5);

    xlim([0.7 1.3]);
    xticks(1);
    xticklabels({summaryAnimalID{s}});

    ylabel('Skew', 'FontSize', 16);

    title(sprintf('%s (%d / %d)', ...
        summaryAnimalID{s}, summaryNSig(s), summaryNTotal(s)), ...
        'FontSize', 16, 'FontWeight', 'bold');

    set(gca, 'FontSize', 15);
    box off;
    grid on;

    if s == 1 && isgraphics(hCI)
        lgd = legend([hDot hCI], ...
            {'Actual Skew', '95% Null Skew Interval'}, ...
            'Location', 'southoutside', 'FontSize', 16);
        lgd.Layout.Tile = 'south';
    end
end

title(t, 'Summary Of Skew Relative To Null Distribution Across Sessions', ...
    'FontSize', 18, 'FontWeight', 'bold');

end

function edges = makeLagEdges(x)
x = x(~isnan(x) & isfinite(x));

if isempty(x)
    edges = -0.201:0.001:0.201;
    return;
end

xmin = floor(min(x)*1000)/1000;
xmax = ceil(max(x)*1000)/1000;

if xmin == xmax
    xmin = xmin - 0.001;
    xmax = xmax + 0.001;
end

edges = xmin:0.001:xmax;

if numel(edges) < 2
    edges = [xmin-0.001 xmax+0.001];
end
end
