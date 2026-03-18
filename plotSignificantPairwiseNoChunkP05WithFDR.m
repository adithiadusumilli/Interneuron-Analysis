function plotSignificantPairwiseNoChunkP05WithFDR(combinedMatFile, alpha, nNullDraws, corrThresh)
% plots significant pairwise NO-CHUNK xcorr results using:
%   p < alpha AND peakCorr > corrThresh
% also reports the storey-estimated FDR at that threshold

% old plotting style kept for:
%   - heatmap
%   - scatterhist
%   - peak lag histogram

% new plot:
%   - null skew distribution vs actual skew

% required variables in combinedMatFile:
%   allSessions.sessions(s).peakCorrMatAll
%   allSessions.sessions(s).peakLagMatAll
%   allSessions.sessions(s).nullCorrMatAllShifts
%   allSessions.sessions(s).nInt
%   allSessions.sessions(s).nPyr
%   allSessions.sessions(s).nAll
%   allSessions.sessions(s).baseDir

% j run:
% plotSignificantPairwiseNoChunkP05WithFDR("C:\Users\mirilab\Documents\GlobusTransfer\pairwise_nochunk_allPairs_ALL_SESSIONS_COMBINED.mat", 0.05, 100, 0.1)

arguments
    combinedMatFile (1,1) string
    alpha (1,1) double = 0.05
    nNullDraws (1,1) double = 100
    corrThresh (1,1) double = 0.1
end

if exist('mafdr', 'file') ~= 2
    error('mafdr is not on the matlab path / bioinformatics toolbox is unavailable.');
end

S = load(combinedMatFile, 'allSessions');
sessions = S.allSessions.sessions;
numSessions = numel(sessions);

summaryActualSkew = nan(numSessions,1);
summaryNullCI = nan(numSessions,2);
summaryFDR = nan(numSessions,1);
summaryAnimalID = cell(1,numSessions);

rng(0);

for sess = 1:numSessions
    peakCorrs = sessions(sess).peakCorrMatAll;
    peakLags  = sessions(sess).peakLagMatAll;
    nullXC    = sessions(sess).nullCorrMatAllShifts;

    nInt = sessions(sess).nInt;
    nPyr = sessions(sess).nPyr;
    nAll = sessions(sess).nAll;

    animalID = regexp(sessions(sess).baseDir, 'D\d+', 'match', 'once');
    if isempty(animalID)
        animalID = sprintf('Session%d', sess);
    end
    summaryAnimalID{sess} = animalID;

    if isempty(peakCorrs) || isempty(peakLags) || isempty(nullXC)
        fprintf('sess %d: empty data, skipping\n', sess);
        continue;
    end

    [rows, cols] = getIntPyrPairs(nInt, nPyr);
    actual = computePairStatsP05(peakCorrs, peakLags, nullXC, rows, cols, alpha, corrThresh);

    fprintf('\n=== %s no-chunk ===\n', animalID);
    fprintf('pairs: %d total | %d significant at p<%.3f & corr>%.3f | est FDR = %.4f\n', ...
        actual.nPairsNominal, actual.nSig, alpha, corrThresh, actual.fdrAtThreshold);

    % -------- heatmap (old style) --------
    sigMat = nan(nInt, nPyr);
    for k = 1:numel(actual.rows)
        r = actual.rows(k);
        c = actual.cols(k) - nInt;
        if actual.sigMask(k)
            sigMat(r,c) = actual.lagVals(k);
        end
    end

    figure('Name', sprintf('%s significant pairwise no-chunk peak lag', animalID), 'Color', 'w');
    hImg = imagesc(sigMat);
    set(hImg, 'AlphaData', ~isnan(sigMat));
    set(gca, 'Color', 'k');
    axis xy;
    colormap(parula);
    cb = colorbar;
    ylabel(cb, 'peak lag (s)');
    xlabel('pyramidal neurons');
    ylabel('interneurons');
    title(sprintf('%s – sig pairs: %d / %d | p < %.3f & corr > %.2f | est FDR = %.3f', ...
        animalID, actual.nSig, actual.nPairsNominal, alpha, corrThresh, actual.fdrAtThreshold));
    set(gca, 'TickDir', 'out');
    box off;

    % -------- lag histogram --------
    figure('Color','w');
    if ~isempty(actual.sigLagVec)
        lagBinEdgesSec = makeLagEdges(actual.sigLagVec);
        histogram(actual.sigLagVec, 'BinEdges', lagBinEdgesSec, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    else
        histogram(nan, 'BinEdges', -0.201:0.001:0.201, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    end
    xlabel('peak lag (s)');
    ylabel('count');
    title(sprintf('%s – significant peak lags (n=%d / %d) | p < %.3f & corr > %.2f | est FDR = %.3f', ...
        animalID, actual.nSig, actual.nPairsNominal, alpha, corrThresh, actual.fdrAtThreshold));
    grid on;

    % -------- scatterhist --------
    if ~isempty(actual.sigLagVec)
        figure('Color','w');
        scatterhist(actual.sigLagVec(:), actual.sigCorrVec(:), 'Direction', 'out', 'Marker', '.');
        xlabel('peak lag (s)');
        ylabel('peak correlation');
        title(sprintf('%s – significant pairs (n=%d / %d) | p < %.3f & corr > %.2f | est FDR = %.3f', ...
            animalID, actual.nSig, actual.nPairsNominal, alpha, corrThresh, actual.fdrAtThreshold));
    end

    % -------- skew null distribution --------
    [poolRows, poolCols] = getAllUpperPairs(nAll);
    nullSkews = nan(nNullDraws,1);

    for r = 1:nNullDraws
        drawIdx = randperm(numel(poolRows), actual.nPairsNominal);
        drawRows = poolRows(drawIdx);
        drawCols = poolCols(drawIdx);

        nullDraw = computePairStatsP05(peakCorrs, peakLags, nullXC, drawRows, drawCols, alpha, corrThresh);
        nullSkews(r) = nullDraw.skew;
    end

    validNullSkews = nullSkews(~isnan(nullSkews) & isfinite(nullSkews));
    if isempty(validNullSkews) || isnan(actual.skew)
        nullCI = [NaN NaN];
        skewP = NaN;
    else
        nullCI = prctile(validNullSkews, [2.5 97.5]);
        skewP = (sum(validNullSkews >= actual.skew) + 1) / (numel(validNullSkews) + 1);
    end

    summaryActualSkew(sess) = actual.skew;
    summaryNullCI(sess,:) = nullCI;
    summaryFDR(sess) = actual.fdrAtThreshold;

    figure('Color','w');
    if ~isempty(validNullSkews)
        histogram(validNullSkews, 30, 'FaceAlpha', 0.8, 'EdgeColor', 'none'); hold on;
        xline(actual.skew, 'r-', 'LineWidth', 2);
        xline(nullCI(1), '--k', 'LineWidth', 1.2);
        xline(nullCI(2), '--k', 'LineWidth', 1.2);
        legend({'Null skews', 'Actual skew', 'Null 95% CI', ''}, 'Location', 'best');
    else
        histogram(nan, 30, 'FaceAlpha', 0.8, 'EdgeColor', 'none'); hold on;
        xline(actual.skew, 'r-', 'LineWidth', 2);
    end
    xlabel('Skew = (mean - median) / std');
    ylabel('Count');
    title(sprintf('%s – skew null | actual = %.4f | p = %.4f | est FDR = %.3f', ...
        animalID, actual.skew, skewP, actual.fdrAtThreshold));
    grid on;
end

% -------- summary plot --------
figure('Name', 'no-chunk skew summary', 'Color', 'w');
tiledlayout(1, numSessions, 'TileSpacing', 'compact', 'Padding', 'compact');

for s = 1:numSessions
    nexttile; hold on
    if ~any(isnan(summaryNullCI(s,:)))
        line([1 1], summaryNullCI(s,:), 'Color', [0.5 0.5 0.5], 'LineWidth', 3);
        plot(1, summaryNullCI(s,1), '_', 'Color', [0.3 0.3 0.3], 'MarkerSize', 18, 'LineWidth', 2);
        plot(1, summaryNullCI(s,2), '_', 'Color', [0.3 0.3 0.3], 'MarkerSize', 18, 'LineWidth', 2);
    end
    plot(1, summaryActualSkew(s), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
    yline(0, 'k:');
    xlim([0.7 1.3]);
    xticks(1);
    xticklabels({summaryAnimalID{s}});
    ylabel('Skew');
    title(sprintf('%s\nFDR@(p<%.2f & corr>%.2f) = %.3f', summaryAnimalID{s}, alpha, corrThresh, summaryFDR(s)));
    box off
    grid on
end

sgtitle('No-chunk summary: actual skew vs null 95% CI');

end

function out = computePairStatsP05(peakCorrMat, peakLagMat, nullCorrMat, rows, cols, alpha, corrThresh)
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

    if isnan(thisReal) || isnan(thisLag) || numel(thisNull) < 10
        pVals(i) = NaN;
        continue
    end

    pVals(i) = (sum(thisNull >= thisReal) + 1) / (numel(thisNull) + 1);
end

validP = ~isnan(pVals) & isfinite(pVals);
sigMask = false(size(pVals));
sigMask(validP) = (pVals(validP) < alpha) & (realVals(validP) > corrThresh);

qVals = nan(size(pVals));
if any(validP)
    [~, qTmp] = mafdr(pVals(validP));
    qVals(validP) = qTmp;
end

if any(sigMask)
    fdrAtThreshold = mean(qVals(sigMask), 'omitnan');
else
    fdrAtThreshold = NaN;
end

sigLagVec = lagVals(sigMask);
sigLagVec = sigLagVec(~isnan(sigLagVec) & isfinite(sigLagVec));

sigCorrVec = realVals(sigMask);
sigCorrVec = sigCorrVec(~isnan(sigCorrVec) & isfinite(sigCorrVec));

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
out.fdrAtThreshold = fdrAtThreshold;
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
