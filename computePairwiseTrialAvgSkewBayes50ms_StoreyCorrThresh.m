function computePairwiseTrialAvgSkewBayes50ms_StoreyCorrThresh(alpha, corrThresh, nNullSkewDraws)
% Pairwise trial-avg chunked Bayes-style skew analysis

% REAL:
%   original int-pyr pairwise peak lags/corrs
%   -> Storey/mafdr + corrThresh significant pairs
%   -> actual skew

% H0 NULL:
%   regular all-vs-all null pair peak lags/corrs
%   -> Storey/mafdr + corrThresh
%   -> null skew distribution

% H50 NULL:
%   shift interneuron trial-avg traces back by 50 ms
%   -> recompute all-vs-all null pair peak lags/corrs
%   -> Storey/mafdr + corrThresh
%   -> shifted null skew distribution

% Evidence ratio: evidenceRatio_H0_over_H50 = pValH0 / pValH50

if nargin < 1 || isempty(alpha), alpha = 0.05; end
if nargin < 2 || isempty(corrThresh), corrThresh = 0.05; end
if nargin < 3 || isempty(nNullSkewDraws), nNullSkewDraws = 1000; end

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

rng(1);

results = struct();

for sessInd = 1:numel(baseDirs)

    baseDir = baseDirs{sessInd};
    outDir = fullfile(baseDir, 'quest_runs');

    fprintf('\n=============================\n');
    fprintf('processing %s\n', baseDir);
    fprintf('=============================\n');

    %% ---------------- build trial-avg traces ----------------
    [allAvg, neuronType, pyrIdx, intIdx, tAxis, binSize, ok] = ...
        buildTrialAvgAllTraces(baseDir, channelsToUse, doBaselineNorm);

    if ~ok
        warning('could not build trial-avg traces for %s', baseDir);
        continue;
    end

    T = numel(tAxis);
    [~, cIdx] = min(abs(tAxis - 0));
    chunkStart = cIdx - chunkHalf;
    chunkEnd   = cIdx + chunkHalf;

    reqMaxLagBins = round(maxLagSecs / binSize);
    maxLagLeft  = chunkStart - 1;
    maxLagRight = T - chunkEnd;
    maxLagBins = min([reqMaxLagBins, maxLagLeft, maxLagRight]);
    lags = -maxLagBins:maxLagBins;

    %% ---------------- real int-pyr pairwise values ----------------
    [realLagVals, realCorrVals] = computeRealIntPyrPairVals( ...
        allAvg, pyrIdx, intIdx, lags, chunkStart, chunkEnd, T, binSize);

    nPairsNominal = numel(realLagVals);

    %% ---------------- regular H0 null from saved null files ----------------
    [nullLagH0_raw, nullCorrH0_raw] = loadRegularNullPairs(outDir, sessInd);

    if isempty(nullLagH0_raw) || isempty(nullCorrH0_raw)
        warning('missing regular null pair values for %s', baseDir);
        continue;
    end

    %% ---------------- H50 null: shift int traces back by 50 ms and recompute all-vs-all ----------------
    lagBin50 = round(lagShiftSec / binSize);

    allAvgShift50 = allAvg;
    for ii = 1:numel(intIdx)
        r = intIdx(ii);
        shiftedTrace = nan(size(allAvg(r,:)));
        shiftedTrace(1:end-lagBin50) = allAvg(r, 1+lagBin50:end);
        allAvgShift50(r,:) = shiftedTrace;
    end

    [nullLagH50_raw, nullCorrH50_raw] = computeAllOrderedPairNullVals( ...
        allAvgShift50, lags, chunkStart, chunkEnd, T, binSize);

    %% ---------------- actual significant skew under original data ----------------
    pReal = corrPvalsFromNull(realCorrVals, nullCorrH0_raw);
    qReal = mafdrSafe(pReal);

    sigRealMask = qReal <= alpha & realCorrVals > corrThresh;

    sigRealLags = realLagVals(sigRealMask);
    actualSkew = computeSkewMetric(sigRealLags);

    %% ---------------- H0 null skew distribution ----------------
    nullSkewH0 = computeNullSkewDistribution( ...
        nullLagH0_raw, nullCorrH0_raw, nullCorrH0_raw, ...
        nPairsNominal, nNullSkewDraws, alpha, corrThresh);

    %% ---------------- H50 null skew distribution ----------------
    nullSkewH50 = computeNullSkewDistribution( ...
        nullLagH50_raw, nullCorrH50_raw, nullCorrH50_raw, ...
        nPairsNominal, nNullSkewDraws, alpha, corrThresh);

    validH0 = nullSkewH0(~isnan(nullSkewH0) & isfinite(nullSkewH0));
    validH50 = nullSkewH50(~isnan(nullSkewH50) & isfinite(nullSkewH50));

    %% ---------------- p-values and evidence ratio ----------------
    pValH0 = (sum(abs(validH0) >= abs(actualSkew)) + 1) / (numel(validH0) + 1);
    pValH50 = (sum(abs(validH50) >= abs(actualSkew)) + 1) / (numel(validH50) + 1);

    evidenceRatio_H0_over_H50 = pValH0 / pValH50;

    %% ---------------- store ----------------
    results(sessInd).baseDir = baseDir;
    results(sessInd).sessInd = sessInd;
    results(sessInd).alpha = alpha;
    results(sessInd).corrThresh = corrThresh;
    results(sessInd).nNullSkewDraws = nNullSkewDraws;
    results(sessInd).lagShiftSec = lagShiftSec;

    results(sessInd).nPairsNominal = nPairsNominal;
    results(sessInd).nSigReal = sum(sigRealMask);
    results(sessInd).sigRealMask = sigRealMask;
    results(sessInd).sigRealLags = sigRealLags;
    results(sessInd).actualSkew = actualSkew;

    results(sessInd).nullSkewH0 = nullSkewH0;
    results(sessInd).nullSkewH50 = nullSkewH50;

    results(sessInd).pValH0 = pValH0;
    results(sessInd).pValH50 = pValH50;
    results(sessInd).evidenceRatio_H0_over_H50 = evidenceRatio_H0_over_H50;

    results(sessInd).realLagVals = realLagVals;
    results(sessInd).realCorrVals = realCorrVals;
    results(sessInd).pReal = pReal;
    results(sessInd).qReal = qReal;

    fprintf('actual skew = %.6f\n', actualSkew);
    fprintf('significant real pairs = %d / %d\n', sum(sigRealMask), nPairsNominal);
    fprintf('pValH0 = %.6f | pValH50 = %.6f\n', pValH0, pValH50);
    fprintf('evidenceRatio H0/H50 = %.6f\n', evidenceRatio_H0_over_H50);
end

outFile = '/home/asa7288/pairwiseTrialAvgSkewBayes50ms_StoreyCorrThresh.mat';
save(outFile, 'results', '-v7.3');

fprintf('\nsaved:\n%s\n', outFile);

%% ---------------- quick plots ----------------
plotSkewBayesResults(results);

end

%% ========================================================================
%% helper functions
%% ========================================================================

function [allAvg, neuronType, pyrIdx, intIdx, tAxis, binSize, ok] = ...
    buildTrialAvgAllTraces(baseDir, channelsToUse, doBaselineNorm)

ok = false;
allAvg = [];
neuronType = [];
pyrIdx = [];
intIdx = [];
tAxis = [];
binSize = 0.001;

S = load(fullfile(baseDir, 'EMG_Neural_AllChannels.mat'), ...
    'pyrCxWinCell','intCxWinCell','tAxis');

tAxis = S.tAxis(:)';

validCh = [];
nPyrList = [];
nIntList = [];

for ci = 1:numel(channelsToUse)
    ch = channelsToUse(ci);
    pyrWin = S.pyrCxWinCell{ch};
    intWin = S.intCxWinCell{ch};

    if isempty(pyrWin) || isempty(intWin)
        continue;
    end

    validCh(end+1) = ch;
    nPyrList(end+1) = size(pyrWin,2);
    nIntList(end+1) = size(intWin,2);
end

if isempty(validCh)
    return;
end

nPyr_ref = min(nPyrList);
nInt_ref = min(nIntList);

pyrAll = [];
intAll = [];

for ci = 1:numel(validCh)
    ch = validCh(ci);

    pyrWin = S.pyrCxWinCell{ch};
    intWin = S.intCxWinCell{ch};

    nEvt = min(size(pyrWin,1), size(intWin,1));
    pyrWin = pyrWin(1:nEvt,1:nPyr_ref,:);
    intWin = intWin(1:nEvt,1:nInt_ref,:);

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

pyrAvg = squeeze(mean(pyrAll, 1, 'omitnan'));
intAvg = squeeze(mean(intAll, 1, 'omitnan'));

allAvg = cat(1, pyrAvg, intAvg);

pyrIdx = 1:nPyr_ref;
intIdx = (nPyr_ref+1):(nPyr_ref+nInt_ref);
neuronType = [zeros(1,nPyr_ref), ones(1,nInt_ref)];

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

function [lagVals, corrVals] = computeRealIntPyrPairVals(allAvg, pyrIdx, intIdx, lags, chunkStart, chunkEnd, T, binSize)

nPairs = numel(pyrIdx) * numel(intIdx);
lagVals = nan(nPairs,1);
corrVals = nan(nPairs,1);

w = 0;

for ii = 1:numel(intIdx)
    intNeuron = intIdx(ii);

    for pp = 1:numel(pyrIdx)
        pyrNeuron = pyrIdx(pp);

        w = w + 1;

        % use pyr as first trace and int as shifted trace
        % positive lag means interneuron leads pyramidal
        [~, peakLagSec, peakCorr] = onePairLagSweep( ...
            allAvg(pyrNeuron,:), allAvg(intNeuron,:), ...
            lags, chunkStart, chunkEnd, T, binSize);

        lagVals(w) = peakLagSec;
        corrVals(w) = peakCorr;
    end
end

valid = ~isnan(lagVals) & ~isnan(corrVals) & isfinite(lagVals) & isfinite(corrVals);
lagVals = lagVals(valid);
corrVals = corrVals(valid);

end

function [nullLag, nullCorr] = loadRegularNullPairs(outDir, sessInd)

nullFiles = dir(fullfile(outDir, ...
    sprintf('pairwiseChunkedTrialAvgXCorr_ALLPAIRS_sess%02d_permNull_*.mat', sessInd)));

nullLag = [];
nullCorr = [];

for k = 1:numel(nullFiles)
    D = load(fullfile(outDir, nullFiles(k).name), ...
        'nullPeakLagVec','nullPeakCorrVec');

    if isfield(D, 'nullPeakLagVec')
        nullLag = [nullLag; D.nullPeakLagVec(:)];
    end

    if isfield(D, 'nullPeakCorrVec')
        nullCorr = [nullCorr; D.nullPeakCorrVec(:)];
    end
end

valid = ~isnan(nullLag) & ~isnan(nullCorr) & isfinite(nullLag) & isfinite(nullCorr);
nullLag = nullLag(valid);
nullCorr = nullCorr(valid);

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
    bEnd = chunkEnd - L;

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

function pVals = corrPvalsFromNull(realCorrVals, nullCorrDist)

nullCorrDist = nullCorrDist(~isnan(nullCorrDist) & isfinite(nullCorrDist));
pVals = nan(size(realCorrVals));

for k = 1:numel(realCorrVals)
    pVals(k) = (sum(nullCorrDist >= realCorrVals(k)) + 1) / (numel(nullCorrDist) + 1);
end

end

function qVals = mafdrSafe(pVals)

pVals = pVals(:);
qVals = nan(size(pVals));

valid = ~isnan(pVals) & isfinite(pVals);

if nnz(valid) < 1
    return;
end

try
    qVals(valid) = mafdr(pVals(valid));
catch
    qVals(valid) = bhFDR(pVals(valid));
end

end

function q = bhFDR(p)

p = p(:);
[ps, sortIdx] = sort(p);
m = numel(p);

qs = ps .* m ./ (1:m)';
qs = flipud(cummin(flipud(qs)));
qs(qs > 1) = 1;

q = nan(size(p));
q(sortIdx) = qs;

end

function nullSkews = computeNullSkewDistribution(nullLagRaw, nullCorrRaw, nullCorrDistForP, ...
    nPairsNominal, nNullSkewDraws, alpha, corrThresh)

nullLagRaw = nullLagRaw(:);
nullCorrRaw = nullCorrRaw(:);

valid = ~isnan(nullLagRaw) & ~isnan(nullCorrRaw) & isfinite(nullLagRaw) & isfinite(nullCorrRaw);
nullLagRaw = nullLagRaw(valid);
nullCorrRaw = nullCorrRaw(valid);

nAvailable = numel(nullLagRaw);
nullSkews = nan(nNullSkewDraws,1);

if nAvailable < 1
    return;
end

for d = 1:nNullSkewDraws

    idx = randi(nAvailable, nPairsNominal, 1);

    lagDraw = nullLagRaw(idx);
    corrDraw = nullCorrRaw(idx);

    pDraw = corrPvalsFromNull(corrDraw, nullCorrDistForP);
    qDraw = mafdrSafe(pDraw);

    sigMask = qDraw <= alpha & corrDraw > corrThresh;

    sigLags = lagDraw(sigMask);

    nullSkews(d) = computeSkewMetric(sigLags);
end

end

function sk = computeSkewMetric(lagVals)

lagVals = lagVals(~isnan(lagVals) & isfinite(lagVals));

if numel(lagVals) < 3
    sk = NaN;
    return;
end

sd = std(lagVals, 'omitnan');

if isnan(sd) || sd == 0
    sk = NaN;
else
    sk = (mean(lagVals, 'omitnan') - median(lagVals, 'omitnan')) / sd;
end

end

function plotSkewBayesResults(results)

nSess = numel(results);

animalIDs = strings(nSess,1);
actualSkew = nan(nSess,1);
pH0 = nan(nSess,1);
pH50 = nan(nSess,1);
BF = nan(nSess,1);

for s = 1:nSess
    [~, tag] = fileparts(results(s).baseDir);
    animalIDs(s) = string(tag);
    actualSkew(s) = results(s).actualSkew;
    pH0(s) = results(s).pValH0;
    pH50(s) = results(s).pValH50;
    BF(s) = results(s).evidenceRatio_H0_over_H50;
end

figure('Color','w','Name','Pairwise Skew Bayes Evidence Ratio');
bar(BF);
hold on;
yline(1,'k--','LineWidth',1.5);
xticks(1:nSess);
xticklabels(animalIDs);
ylabel('Evidence Ratio H0 / H50');
xlabel('Animal');
title('Pairwise Skew Bayes Evidence Ratio');
box off;
set(gca,'FontSize',14,'LineWidth',1,'TickDir','out');

figure('Color','w','Name','Pairwise Skew Bayes p-values');
bar([pH0 pH50]);
xticks(1:nSess);
xticklabels(animalIDs);
ylabel('p-value');
xlabel('Animal');
legend({'H0 regular null skew','H50 shifted-int null skew'}, 'Location','best');
title('Pairwise Skew p-values');
box off;
set(gca,'FontSize',14,'LineWidth',1,'TickDir','out');

figure('Color','w','Name','Actual Pairwise Skew');
bar(actualSkew);
hold on;
yline(0,'k--','LineWidth',1.5);
xticks(1:nSess);
xticklabels(animalIDs);
ylabel('Actual skew');
xlabel('Animal');
title('Actual Skew From Significant Int-Pyr Pairs');
box off;
set(gca,'FontSize',14,'LineWidth',1,'TickDir','out');

end
