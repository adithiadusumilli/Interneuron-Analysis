function plotSavedIntPyrLagSkewMAFDR_Chunked(resultsFile, combinedMatFile)
% plots saved chunked pairwise significance results from
% runIntPyrLagSkewMAFDRPermutationTest

% expected inputs:
%   resultsFile = output .mat from runIntPyrLagSkewMAFDRPermutationTest for analysisType = "chunked"
%   combinedMatFile = original combined chunked pairwise file containing all_nullCorrMat_allShifts

% figures per session:
%   1. significant pairwise peak lag map
%   2. distribution of significant peak lags
%   3. peak correlation relative to shift control
%   4. peak lag versus peak correlation
%   5. skew relative to null distribution

% summary figure -- actual skew versus null 95% interval across sessions

% j run: plotSavedIntPyrLagSkewMAFDR_Chunked("C:\Users\mirilab\Documents\GlobusTransfer\pairwiseChunked_ALLPAIRS_ALLSESS_REAL_AND_SHIFTS_COMBINED_intPyrSkewMAFDR_chunked.mat", "C:\Users\mirilab\Documents\GlobusTransfer\pairwiseChunked_ALLPAIRS_ALLSESS_REAL_AND_SHIFTS_COMBINED.mat")

arguments
    resultsFile (1,1) string
    combinedMatFile (1,1) string
end

S = load(resultsFile, 'results');
R = S.results;

C = load(combinedMatFile, 'all_nullCorrMat_allShifts', 'all_xcMat_all', 'all_lags', 'all_binSize');

rawNullShifts = C.all_nullCorrMat_allShifts;

if ~isfield(R, 'sessions') || isempty(R.sessions)
    error('results file does not contain a valid results.sessions field.');
end

if isfield(R, 'analysisType')
    if ~strcmpi(string(R.analysisType), "chunked")
        error('this plotting function is for chunked saved outputs only.');
    end
end

nSess = numel(R.sessions);

if numel(rawNullShifts) ~= nSess
    error('results file and combined chunked file do not have the same number of sessions.');
end

summaryActualSkew = nan(nSess,1);
summaryNullCI = nan(nSess,2);
summaryAnimalID = cell(1,nSess);
summaryNSig = nan(nSess,1);
summaryNTotal = nan(nSess,1);

nullColor = [0 0.4470 0.7410];

for s = 1:nSess
    sess = R.sessions{s};
    rawNullXC = rawNullShifts{s};

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
    box off;

    % ================= EXAMPLE PAIR CORR VS LAG RELATIVE TO SHIFT CONTROL =================
    figure('Name', sprintf('%s Example Pair Corr Versus Lag', animalID), 'Color', 'w');
    hold on;

    hCurve = gobjects(1);
    hPeakLag = gobjects(1);
    hPeak = gobjects(1);
    hCI = gobjects(1);

    sigPairsIdx = find(actual.sigFDRMask(:));

    if ~isempty(sigPairsIdx)
        % choose one significant pair with the largest real peak correlation
        [~, bestLocal] = max(actual.realVals(sigPairsIdx));
        exampleIdx = sigPairsIdx(bestLocal);

        r = actual.rows(exampleIdx);
        c = actual.cols(exampleIdx);

        % get real xc curve and lag vector
        [lagsVec, xcVec] = getChunkedExampleXC(C, s, r, c);

        % get shift-control distribution for this same pair
        thisNull = squeeze(rawNullXC(r,c,:));
        thisNull = thisNull(~isnan(thisNull) & isfinite(thisNull));

        if ~isempty(lagsVec) && ~isempty(xcVec)
            lagsVec = lagsVec(:);
            xcVec = xcVec(:);

            % plot real xc curve
            hCurve = plot(lagsVec, xcVec, '-', ...
                'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);

            % compute peak directly from the real xc curve
            [peakCorrVal, peakIdx] = max(xcVec);
            peakLagVal = lagsVec(peakIdx);

            % mark actual peak lag and peak corr
            hPeakLag = xline(peakLagVal, 'r:', 'LineWidth', 1.8);
            hPeak = plot(peakLagVal, peakCorrVal, 'o', ...
                'MarkerFaceColor', [0.85 0.33 0.10], ...
                'MarkerEdgeColor', [0.85 0.33 0.10], ...
                'MarkerSize', 8);

            % add shift-control bounds
            if ~isempty(thisNull)
                corrCI = prctile(thisNull, [2.5 97.5]);
                yline(corrCI(1), 'k--', 'LineWidth', 1.5);
                hCI = yline(corrCI(2), 'k--', 'LineWidth', 1.5);
            end
        end
    end

    xlabel('Lag (Seconds)', 'FontSize', 14);
    ylabel('Correlation', 'FontSize', 14);

    title(sprintf('%s Example Pair Corr Versus Lag (%d / %d Significant Pairs)', ...
        animalID, actual.nSigFDR, actual.nPairsNominal), ...
        'FontSize', 16, 'FontWeight', 'bold');

    legendHandles = gobjects(0);
    legendLabels = {};

    if isgraphics(hCurve)
        legendHandles(end+1) = hCurve;
        legendLabels{end+1} = 'Example Pair XC';
    end
    if isgraphics(hPeakLag)
        legendHandles(end+1) = hPeakLag;
        legendLabels{end+1} = 'Actual Peak Lag';
    end
    if isgraphics(hPeak)
        legendHandles(end+1) = hPeak;
        legendLabels{end+1} = 'Peak Correlation';
    end
    if isgraphics(hCI)
        legendHandles(end+1) = hCI;
        legendLabels{end+1} = '95% Shift-Control Bounds';
    end
    if ~isempty(legendHandles)
        legend(legendHandles, legendLabels, 'FontSize', 12, 'Location', 'best');
    end

    set(gca, 'FontSize', 13);
    box off;
    
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
        histogram(validNullSkews, 30, 'EdgeColor', 'none', 'FaceColor', nullColor);
        hold on;
        hActual = xline(actual.skew, 'r', 'LineWidth', 2);
        hCI = xline(nullCI(1), 'k--', 'LineWidth', 1.5);
        xline(nullCI(2), 'k--', 'LineWidth', 1.5);

        hNullDummy = patch(nan, nan, nullColor);
        legend([hNullDummy hActual hCI], ...
            {'Null Skew Distribution', 'Actual Skew', '95% Null Skew Interval'}, ...
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
    box off;
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

title(t, 'Summary Of Skew Relative To Null Distribution', ...
    'FontSize', 18, 'FontWeight', 'bold');

end

function [lagsVec, xcVec] = getChunkedExampleXC(C, s, r, c)
lagsVec = [];
xcVec = [];

if isfield(C, 'all_xcMat_all') && numel(C.all_xcMat_all) >= s && ~isempty(C.all_xcMat_all{s})
    xcCube = C.all_xcMat_all{s};
else
    return;
end

xcVec = squeeze(xcCube(r,c,:));

if isfield(C, 'all_lags') && numel(C.all_lags) >= s && ~isempty(C.all_lags{s})
    lagsVec = C.all_lags{s}(:);
elseif isfield(C, 'all_binSize') && numel(C.all_binSize) >= s && ~isempty(C.all_binSize{s})
    binSize = C.all_binSize{s};
    nL = numel(xcVec);
    mid = ceil(nL/2);
    lagsVec = ((1:nL) - mid) * binSize;
else
    lagsVec = [];
end
end

function out = computePairStatsQCorr(peakCorrMat, peakLagMat, nullCorrMat, rows, cols, alpha, corrThresh)
nPairs = numel(rows);
realVals = nan(nPairs,1);
lagVals = nan(nPairs,1);
pVals = nan(nPairs,1);

for i = 1:nPairs
    r = rows(i);
    c = cols(i);

    thisReal = peakCorrMat(r,c);
    thisLag = peakLagMat(r,c);
    thisNull = squeeze(nullCorrMat(r,c,:));
    thisNull = thisNull(~isnan(thisNull) & isfinite(thisNull));

    realVals(i) = thisReal;
    lagVals(i) = thisLag;

    if isnan(thisReal) || ~isfinite(thisReal) || isnan(thisLag) || ~isfinite(thisLag) || numel(thisNull) < 10
        pVals(i) = NaN;
        continue
    end

    pVals(i) = (sum(thisNull >= thisReal) + 1) / (numel(thisNull) + 1);
end

validP = ~isnan(pVals) & isfinite(pVals);

qVals = nan(size(pVals));
if any(validP)
    [~, qTmp] = mafdr(pVals(validP));
    qVals(validP) = qTmp;
end

sigMask = false(size(pVals));
if any(validP)
    sigMask(validP) = (qVals(validP) <= alpha) & (realVals(validP) > corrThresh);
end

sigLagVec = lagVals(sigMask);
sigLagVec = sigLagVec(~isnan(sigLagVec) & isfinite(sigLagVec));

sigCorrVec = realVals(sigMask);
sigCorrVec = sigCorrVec(~isnan(sigCorrVec) & isfinite(sigCorrVec));

if any(sigMask)
    meanQSelected = mean(qVals(sigMask), 'omitnan');
else
    meanQSelected = NaN;
end

out.rows = rows;
out.cols = cols;
out.realVals = realVals;
out.lagVals = lagVals;
out.pVals = pVals;
out.qVals = qVals;
out.sigMask = sigMask;
out.nPairsNominal = nPairs;
out.nValidTests = nnz(validP);
out.nSig = nnz(sigMask);
out.meanQSelected = meanQSelected;
out.sigLagVec = sigLagVec;
out.sigCorrVec = sigCorrVec;
out.skew = computeSkew(sigLagVec);
end

function [rows, cols] = getIntPyrPairs(nInt, nPyr)
[rr, cc] = ndgrid(1:nInt, nInt + (1:nPyr));
rows = rr(:);
cols = cc(:);
end

function [rows, cols] = getAllUpperPairs(nAll)
mask = triu(true(nAll), 1);
[rows, cols] = find(mask);
end

function skewVal = computeSkew(x)
x = x(~isnan(x) & isfinite(x));
if numel(x) < 2
    skewVal = NaN;
    return;
end
sd = std(x, 0, 'omitnan');
if sd == 0 || isnan(sd)
    skewVal = NaN;
    return;
end
skewVal = (mean(x, 'omitnan') - median(x, 'omitnan')) / sd;
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
