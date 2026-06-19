function computePairwiseTrialAvgSkewBayes50ms_StoreyCorrThresh(alpha, corrThresh, nNullDraws)
% pairwise trial-averaged chunked Bayes using the skew metric

% goal:
%   1) Get real skew from significant int-pyr pairs
%   2) Build regular H0 null skew dist
%   3) Build H50 null skew dist after shifting interneuron traces by 50 ms
%   4) Compare same real skew to both null skew dists

% Important:
%   - keeps analysis trial-averaged
%   - Storey/mafdr + corrThresh is applied before calculating skew
%   - 50 ms shift happens before recomputing pairwise peak lags/corrs

if nargin < 1 || isempty(alpha), alpha = 0.05; end
if nargin < 2 || isempty(corrThresh), corrThresh = 0.05; end
if nargin < 3 || isempty(nNullDraws), nNullDraws = 100; end

baseDirs = {
    '/home/asa7288/Transfer/D026', ...
    '/home/asa7288/Transfer/D020', ...
    '/home/asa7288/Transfer/D024', ...
    '/home/asa7288/Transfer/D043'
};

channelsToUse  = 1:4;
chunkHalf = 300;
maxLagSecs = 0.2;
doBaselineNorm = true;
lagShiftSec = 0.050;

rng(0);

results = struct();
results.alpha = alpha;
results.corrThresh = corrThresh;
results.nNullDraws = nNullDraws;
results.lagShiftSec = lagShiftSec;
results.analysisNote = 'trial-avg pairwise skew Bayes; significant pairs selected by Storey/mafdr + corrThresh';

for sessInd = 1:numel(baseDirs)

    baseDir = baseDirs{sessInd};
    animalID = regexp(baseDir, 'D\d+', 'match', 'once');

    fprintf('\n=============================\n');
    fprintf('processing %s\n', animalID);
    fprintf('=============================\n');

    %% build the trial-averaged traces I actually want to analyze
    [allAvg, nInt, nPyr, nAll, tAxis, binSize, ok] = ...
        buildTrialAvgAllTraces(baseDir, channelsToUse, doBaselineNorm);

    if ~ok
        warning('could not build trial-avg traces for %s', baseDir);
        continue;
    end

    T = numel(tAxis);
    [~, cIdx] = min(abs(tAxis - 0));

    chunkStart = cIdx - chunkHalf;
    chunkEnd = cIdx + chunkHalf;

    reqMaxLagBins = round(maxLagSecs / binSize);
    maxLagLeft = chunkStart - 1;
    maxLagRight = T - chunkEnd;
    maxLagBins = min([reqMaxLagBins, maxLagLeft, maxLagRight]);

    lags = -maxLagBins:maxLagBins;

    %% H0/original pairwise matrices from trial-averaged traces
    [peakCorrH0, peakLagH0] = computeAllPairMatrices(allAvg, lags, chunkStart, chunkEnd, T, binSize);

    %% H50 pairwise matrices: shift interneuron traces by 50 ms first, then recompute pairwise
    allAvgH50 = shiftInterneuronsBy50ms(allAvg, nPyr, lagShiftSec, binSize);

    [peakCorrH50, peakLagH50] = computeAllPairMatrices(allAvgH50, lags, chunkStart, chunkEnd, T, binSize);

    %% load regular null correlation values from the saved pairwise trial-avg null jobs
    [nullCorrH0, okNull] = loadTrialAvgNullCorrValues(baseDir, sessInd);

    if ~okNull
        warning('missing null corr files for %s', animalID);
        continue;
    end

    %% recompute H50 null correlations from the shifted-int trial-avg traces
    [~, nullCorrH50] = computeAllOrderedPairNullVals(allAvgH50, lags, chunkStart, chunkEnd, T, binSize);

    %% real int-pyr pair set
    [actualRows, actualCols] = getIntPyrPairs(nInt, nPyr);

    actual = computePairSetStatsMAFDRSkew_GlobalNull(peakCorrH0, peakLagH0, nullCorrH0, actualRows, actualCols, alpha, corrThresh);

    fprintf('%s real: sig pairs = %d / %d | skew = %.6f\n', animalID, actual.nSigFDR, actual.nPairsNominal, actual.skew);

    %% null skew H0: sample random all-vs-all pairs from original trial-avg matrices
    nullSkewH0 = computeNullSkewsFromPairPool(peakCorrH0, peakLagH0, nullCorrH0, nInt, nAll, numel(actualRows), nNullDraws, alpha, corrThresh);

    %% null skew H50: sample random all-vs-all pairs from shifted-int matrices
    nullSkewH50 = computeNullSkewsFromPairPool(peakCorrH50, peakLagH50, nullCorrH50, nInt, nAll, numel(actualRows), nNullDraws, alpha, corrThresh);

    validH0 = nullSkewH0(~isnan(nullSkewH0) & isfinite(nullSkewH0));
    validH50 = nullSkewH50(~isnan(nullSkewH50) & isfinite(nullSkewH50));

    %% Bayes-style p-values from skew distributions
    pValH0 = (sum(abs(validH0) >= abs(actual.skew)) + 1) / (numel(validH0) + 1);
    pValH50 = (sum(abs(validH50) >= abs(actual.skew)) + 1) / (numel(validH50) + 1);

    evidenceRatio_H0_over_H50 = pValH0 / pValH50;

    %% save everything useful
    R = struct();
    R.animalID = animalID;
    R.baseDir = baseDir;
    R.sessInd = sessInd;

    R.nInt = nInt;
    R.nPyr = nPyr;
    R.nAll = nAll;

    R.alpha = alpha;
    R.corrThresh = corrThresh;
    R.nNullDraws = nNullDraws;
    R.lagShiftSec = lagShiftSec;

    R.actual = actual;
    R.nullSkewH0 = nullSkewH0;
    R.nullSkewH50 = nullSkewH50;
    R.validNullSkewH0 = validH0;
    R.validNullSkewH50 = validH50;

    R.pValH0 = pValH0;
    R.pValH50 = pValH50;
    R.evidenceRatio_H0_over_H50 = evidenceRatio_H0_over_H50;

    R.peakCorrH0 = peakCorrH0;
    R.peakLagH0 = peakLagH0;
    R.peakCorrH50 = peakCorrH50;
    R.peakLagH50 = peakLagH50;

    results.sessions{sessInd} = R;

    fprintf('%s Bayes skew: pH0 = %.6f | pH50 = %.6f | ratio H0/H50 = %.6f\n', animalID, pValH0, pValH50, evidenceRatio_H0_over_H50);
end

outFile = '/home/asa7288/pairwiseTrialAvgSkewBayes50ms_StoreyCorrThresh.mat';
save(outFile, 'results', '-v7.3');

fprintf('\nsaved:\n%s\n', outFile);

end

%% ========================================================================
%% helpers
%% ========================================================================

function [allAvg, nInt, nPyr, nAll, tAxis, binSize, ok] = buildTrialAvgAllTraces(baseDir, channelsToUse, doBaselineNorm)

ok = false;
binSize = 0.001;

S = load(fullfile(baseDir, 'EMG_Neural_AllChannels.mat'), ...
    'pyrCxWinCell','intCxWinCell','tAxis');

tAxis = S.tAxis(:)';

validCh = [];
nPyrList = [];
nIntList = [];

for ci = 1:numel(channelsToUse)
    ch = channelsToUse(ci);

    if isempty(S.pyrCxWinCell{ch}) || isempty(S.intCxWinCell{ch})
        continue;
    end

    validCh(end+1) = ch;
    nPyrList(end+1) = size(S.pyrCxWinCell{ch}, 2);
    nIntList(end+1) = size(S.intCxWinCell{ch}, 2);
end

if isempty(validCh)
    allAvg = [];
    nInt = NaN;
    nPyr = NaN;
    nAll = NaN;
    return;
end

nPyr = min(nPyrList);
nInt = min(nIntList);

pyrAll = [];
intAll = [];

for ci = 1:numel(validCh)
    ch = validCh(ci);

    pyrWin = S.pyrCxWinCell{ch};
    intWin = S.intCxWinCell{ch};

    nEvt = min(size(pyrWin,1), size(intWin,1));

    pyrWin = pyrWin(1:nEvt, 1:nPyr, :);
    intWin = intWin(1:nEvt, 1:nInt, :);

    if doBaselineNorm
        pyrWin = subtractTrialBaseline3d(pyrWin, tAxis, -500, -450);
        intWin = subtractTrialBaseline3d(intWin, tAxis, -500, -450);
    end

    pyrAll = cat(1, pyrAll, pyrWin);
    intAll = cat(1, intAll, intWin);
end

nEvt = min(size(pyrAll,1), size(intAll,1));
pyrAll = pyrAll(1:nEvt,:,:);
intAll = intAll(1:nEvt,:,:);

% trial average happens here: events/windows are averaged before pairwise xcorr
pyrAvg = squeeze(mean(pyrAll, 1, 'omitnan'));
intAvg = squeeze(mean(intAll, 1, 'omitnan'));

% ordering is [interneurons, pyramidal] to match my older pairwise code
allAvg = cat(1, intAvg, pyrAvg);

nAll = nInt + nPyr;
ok = true;

end

function Wout = subtractTrialBaseline3d(Win, tAxis, tStart, tEnd)

[~, iStart] = min(abs(tAxis - tStart));
[~, iEnd] = min(abs(tAxis - tEnd));

if iEnd < iStart
    tmp = iStart;
    iStart = iEnd;
    iEnd = tmp;
end

base = mean(Win(:,:,iStart:iEnd), 3, 'omitnan');
Wout = Win - base;

end

function allAvgH50 = shiftInterneuronsBy50ms(allAvg, nPyr, lagShiftSec, binSize)

% allAvg is ordered [int, pyr], so interneurons are the first nInt rows.
% I shift interneuron traces earlier by 50 ms so this represents the H50 model.

lagBin50 = round(lagShiftSec / binSize);

nAll = size(allAvg,1);
nInt = nAll - nPyr;

allAvgH50 = allAvg;

for r = 1:nInt
    shiftedTrace = nan(size(allAvg(r,:)));
    shiftedTrace(1:end-lagBin50) = allAvg(r, 1+lagBin50:end);
    allAvgH50(r,:) = shiftedTrace;
end

end

function [peakCorrMat, peakLagMat] = computeAllPairMatrices(allAvg, lags, chunkStart, chunkEnd, T, binSize)

nAll = size(allAvg,1);

peakCorrMat = nan(nAll, nAll);
peakLagMat = nan(nAll, nAll);

for i = 1:nAll
    for j = (i+1):nAll

        [~, peakLagSec, peakCorr] = onePairLagSweep( ...
            allAvg(i,:), allAvg(j,:), ...
            lags, chunkStart, chunkEnd, T, binSize);

        peakCorrMat(i,j) = peakCorr;
        peakLagMat(i,j) = peakLagSec;
    end
end

end

function [nullLag, nullCorr] = computeAllOrderedPairNullVals(allAvg, lags, chunkStart, chunkEnd, T, binSize)

nAll = size(allAvg,1);
nPairs = nAll * (nAll - 1);

nullLag = nan(nPairs,1);
nullCorr = nan(nPairs,1);

w = 0;

for i = 1:nAll
    for j = 1:nAll

        if i == j
            continue;
        end

        w = w + 1;

        [~, peakLagSec, peakCorr] = onePairLagSweep( ...
            allAvg(i,:), allAvg(j,:), ...
            lags, chunkStart, chunkEnd, T, binSize);

        nullLag(w) = peakLagSec;
        nullCorr(w) = peakCorr;
    end
end

valid = ~isnan(nullLag) & ~isnan(nullCorr) & isfinite(nullLag) & isfinite(nullCorr);
nullLag = nullLag(valid);
nullCorr = nullCorr(valid);

end

function [xc, peakLagSec, peakCorr] = onePairLagSweep(aTrace, bTrace, lags, chunkStart, chunkEnd, T, binSize)

xc = nan(1, numel(lags));

for iL = 1:numel(lags)

    L = lags(iL);

    bStart = chunkStart - L;
    bEnd   = chunkEnd   - L;

    if bStart < 1 || bEnd > T
        continue;
    end

    aseg = aTrace(chunkStart:chunkEnd);
    bseg = bTrace(bStart:bEnd);

    valid = ~isnan(aseg) & ~isnan(bseg);

    if nnz(valid) > 10
        xc(iL) = corr(aseg(valid)', bseg(valid)');
    end
end

[peakCorr, peakIdx] = max(xc);
peakLagSec = lags(peakIdx) * binSize;

end

function [rows, cols] = getIntPyrPairs(nInt, nPyr)

% allAvg is [int, pyr], so int rows are 1:nInt and pyr rows are nInt+1:nAll

[rr, cc] = ndgrid(1:nInt, nInt + (1:nPyr));
rows = rr(:);
cols = cc(:);

end

function [rows, cols] = getAllUpperPairs(nAll)

mask = triu(true(nAll), 1);
[rows, cols] = find(mask);

end

function [nullCorr, ok] = loadTrialAvgNullCorrValues(baseDir, sessInd)

outDir = fullfile(baseDir, 'quest_runs');

nullFiles = dir(fullfile(outDir, ...
    sprintf('pairwiseChunkedTrialAvgXCorr_ALLPAIRS_sess%02d_permNull_*.mat', sessInd)));

nullCorr = [];

for k = 1:numel(nullFiles)
    D = load(fullfile(outDir, nullFiles(k).name), 'nullPeakCorrVec');

    if isfield(D, 'nullPeakCorrVec')
        nullCorr = [nullCorr; D.nullPeakCorrVec(:)];
    end
end

nullCorr = nullCorr(~isnan(nullCorr) & isfinite(nullCorr));
ok = ~isempty(nullCorr);

end

function out = computePairSetStatsMAFDRSkew_GlobalNull(realMat, lagMat, nullCorr, rows, cols, alpha, corrThresh, lagSignVec)

% This mirrors my old Storey/corrThresh function, except here the trial-avg
% pairwise null is a global all-vs-all null correlation distribution.

if nargin < 8 || isempty(lagSignVec)
    lagSignVec = ones(numel(rows),1);
end

nPairs = numel(rows);

realVals = nan(nPairs,1);
lagVals = nan(nPairs,1);
pVals = nan(nPairs,1);

nullCorr = nullCorr(~isnan(nullCorr) & isfinite(nullCorr));

for i = 1:nPairs

    r = rows(i);
    c = cols(i);

    thisReal = realMat(r,c);
    thisLag = lagMat(r,c);

    realVals(i) = thisReal;
    lagVals(i) = lagSignVec(i) * thisLag;

    if isnan(thisReal) || isnan(thisLag) || isempty(nullCorr)
        continue;
    end

    pVals(i) = (sum(nullCorr >= thisReal) + 1) / (numel(nullCorr) + 1);
end

validP = ~isnan(pVals) & isfinite(pVals);

sigUncorr = false(size(pVals));
sigUncorr(validP) = pVals(validP) <= alpha;

qVals = nan(size(pVals));
pFDRVals = nan(size(pVals));

if any(validP)
    [pFDRtmp, qTmp] = mafdr(pVals(validP));
    pFDRVals(validP) = pFDRtmp;
    qVals(validP) = qTmp;
end

sigFDR = false(size(pVals));
sigFDR(validP) = (qVals(validP) <= alpha) & (realVals(validP) > corrThresh);

sigLagVec = lagVals(sigFDR);
sigLagVec = sigLagVec(~isnan(sigLagVec) & isfinite(sigLagVec));

out = struct();
out.rows = rows;
out.cols = cols;
out.realVals = realVals;
out.lagVals = lagVals;
out.pVals = pVals;
out.pFDRVals = pFDRVals;
out.qVals = qVals;
out.qAlpha = alpha;
out.corrThresh = corrThresh;
out.nPairsNominal = nPairs;
out.nValidTests = nnz(validP);
out.nSigUncorr = nnz(sigUncorr);
out.nSigFDR = nnz(sigFDR);
out.sigLagVec = sigLagVec;
out.skew = computeSkew(sigLagVec);
out.sigUncorrMask = sigUncorr;
out.sigFDRMask = sigFDR;
out.lagSignVec = lagSignVec;

end

function nullSkews = computeNullSkewsFromPairPool(peakCorr, peakLag, nullCorr, nInt, nAll, nSamplePairs, nNullDraws, alpha, corrThresh)

[poolRows, poolCols] = getAllUpperPairs(nAll);

if nSamplePairs > numel(poolRows)
    error('number of requested sampled pairs exceeds all upper-triangle pairs');
end

nullSkews = nan(nNullDraws,1);

for d = 1:nNullDraws

    drawIdx = randperm(numel(poolRows), nSamplePairs);

    drawRows = poolRows(drawIdx);
    drawCols = poolCols(drawIdx);

    % I randomly flip the lag sign like in the older code so the sampled
    % all-vs-all null can represent both pair orientations.
    lagSignVec = ones(nSamplePairs,1);
    lagSignVec(rand(nSamplePairs,1) > 0.5) = -1;

    nullDraw = computePairSetStatsMAFDRSkew_GlobalNull( ...
        peakCorr, peakLag, nullCorr, ...
        drawRows, drawCols, alpha, corrThresh, lagSignVec);

    nullSkews(d) = nullDraw.skew;
end

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
else
    skewVal = (mean(x, 'omitnan') - median(x, 'omitnan')) / sd;
end

end
