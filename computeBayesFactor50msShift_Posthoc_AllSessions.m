function computeBayesFactor50msShift_Posthoc_AllSessions()
% computes BF-style evidence ratio:
% BF = pValOriginal / pValShift50
%
% Model 1: original int/pyr timing
% Model 2: interneuron trace shifted back by 50 ms, then xcorr recomputed
%
% run after real and permutation files already exist

baseDirs = {
    'C:\Users\mirilab\Documents\GlobusTransfer\D026', ...
    'C:\Users\mirilab\Documents\GlobusTransfer\D020', ...
    'C:\Users\mirilab\Documents\GlobusTransfer\D024', ...
    'C:\Users\mirilab\Documents\GlobusTransfer\D043'
};

for iDir = 1:numel(baseDirs)

    baseDir = baseDirs{iDir};
    outDir = fullfile(baseDir, 'quest_runs');

    realFile = fullfile(outDir, 'concatCrossCorr_trialavg_chunked_popavg_unperm.mat');

    if ~isfile(realFile)
        warning('missing real file: %s', realFile);
        continue;
    end

    R = load(realFile, ...
        'xc','lags','peakLagSec','pyrAvgTrace','intAvgTrace','binSize','chunkHalf');

    S2 = load(fullfile(baseDir, 'EMG_Neural_AllChannels.mat'), 'tAxis');
    tAxis = S2.tAxis(:)';
    T = numel(tAxis);

    [~, cIdx] = min(abs(tAxis - 0));
    pyrStart = cIdx - R.chunkHalf;
    pyrEnd   = cIdx + R.chunkHalf;

    [bf_H0_over_H1, pValOriginal, pValShift50, ...
        realPeakCorrOriginal, realPeakCorrShift50, ...
        peakLagSecOriginal, peakLagSecShift50, ...
        expectedPeakLagShift50, sanityLagDiff] = computeBFShift50ms( ...
            R.pyrAvgTrace, R.intAvgTrace, R.lags, ...
            pyrStart, pyrEnd, T, R.binSize, R.xc, R.peakLagSec, outDir);

    save(realFile, ...
        'bf_H0_over_H1', ...
        'pValOriginal','pValShift50', ...
        'realPeakCorrOriginal','realPeakCorrShift50', ...
        'peakLagSecOriginal','peakLagSecShift50', ...
        'expectedPeakLagShift50','sanityLagDiff', ...
        '-append');

    fprintf('\n%s\n', baseDir);
    fprintf('BF H0/H1 = %.6f\n', bf_H0_over_H1);
    fprintf('pValOriginal = %.6f | pValShift50 = %.6f\n', pValOriginal, pValShift50);
end

end

function [bf_H0_over_H1, pValOriginal, pValShift50, ...
    realPeakCorrOriginal, realPeakCorrShift50, ...
    peakLagSecOriginal, peakLagSecShift50, ...
    expectedPeakLagShift50, sanityLagDiff] = computeBFShift50ms( ...
    pyrAvgTrace, intAvgTrace, lags, pyrStart, pyrEnd, T, binSize, xcOriginal, peakLagSecOriginal, outDir)

lagBin50ms = round(0.050 / binSize);

% shift interneuron trace back by 50 ms
% with your sign convention, this should shift the peak lag by about +50 ms
intAvgTraceShift50 = nan(size(intAvgTrace));
intAvgTraceShift50(1:end-lagBin50ms) = intAvgTrace(1+lagBin50ms:end);

% recompute xcorr after shifting interneuron trace
[xcShift50, peakLagSecShift50] = lagSweepTrialAverage( ...
    pyrAvgTrace, intAvgTraceShift50, lags, pyrStart, pyrEnd, T, binSize);

realPeakCorrOriginal = max(xcOriginal, [], 'omitnan');
realPeakCorrShift50  = max(xcShift50, [], 'omitnan');

expectedPeakLagShift50 = peakLagSecOriginal + 0.050;
sanityLagDiff = peakLagSecShift50 - expectedPeakLagShift50;

fprintf('\noriginal peak lag = %.6f s | original peak corr = %.6f\n', ...
    peakLagSecOriginal, realPeakCorrOriginal);
fprintf('shifted-int peak lag = %.6f s | shifted-int peak corr = %.6f\n', ...
    peakLagSecShift50, realPeakCorrShift50);
fprintf('expected shifted peak lag approx = %.6f s | diff = %.6f s\n', ...
    expectedPeakLagShift50, sanityLagDiff);

permFiles = dir(fullfile(outDir, 'concatCrossCorr_trialavg_chunked_popavg_perm_*.mat'));

if isempty(permFiles)
    warning('no permutation files found in %s', outDir);
    bf_H0_over_H1 = NaN;
    pValOriginal = NaN;
    pValShift50 = NaN;
    return;
end

permPeakCorrOriginal = nan(numel(permFiles),1);
permPeakCorrShift50  = nan(numel(permFiles),1);

for k = 1:numel(permFiles)

    P = load(fullfile(outDir, permFiles(k).name), ...
        'xc','lags','pyrAvgTrace','intAvgTrace','binSize');

    if isfield(P, 'xc')
        permPeakCorrOriginal(k) = max(P.xc, [], 'omitnan');
    end

    if isfield(P, 'pyrAvgTrace') && isfield(P, 'intAvgTrace')
        intPermShift50 = nan(size(P.intAvgTrace));
        intPermShift50(1:end-lagBin50ms) = P.intAvgTrace(1+lagBin50ms:end);

        [xcPermShift50, ~] = lagSweepTrialAverage( ...
            P.pyrAvgTrace, intPermShift50, lags, pyrStart, pyrEnd, T, binSize);

        permPeakCorrShift50(k) = max(xcPermShift50, [], 'omitnan');
    end
end

permPeakCorrOriginal = permPeakCorrOriginal(~isnan(permPeakCorrOriginal) & isfinite(permPeakCorrOriginal));
permPeakCorrShift50  = permPeakCorrShift50(~isnan(permPeakCorrShift50) & isfinite(permPeakCorrShift50));

% add-one correction
pValOriginal = (sum(permPeakCorrOriginal >= realPeakCorrOriginal) + 1) / ...
               (numel(permPeakCorrOriginal) + 1);

pValShift50 = (sum(permPeakCorrShift50 >= realPeakCorrShift50) + 1) / ...
              (numel(permPeakCorrShift50) + 1);

bf_H0_over_H1 = pValOriginal / pValShift50;

fprintf('pValOriginal = %.6f\n', pValOriginal);
fprintf('pValShift50  = %.6f\n', pValShift50);
fprintf('BF-style ratio H0/H1 = %.6f\n', bf_H0_over_H1);

end

function [xc, peakLagSec] = lagSweepTrialAverage(pyrAvgTrace, intAvgTrace, lags, pyrStart, pyrEnd, T, binSize)

xc = nan(1, numel(lags));

for iL = 1:numel(lags)
    L = lags(iL);

    intStart = pyrStart - L;
    intEnd   = pyrEnd   - L;

    if intStart < 1 || intEnd > T
        continue;
    end

    pseg = pyrAvgTrace(pyrStart:pyrEnd);
    iseg = intAvgTrace(intStart:intEnd);

    valid = ~isnan(pseg) & ~isnan(iseg);

    if nnz(valid) > 10
        xc(iL) = corr(pseg(valid)', iseg(valid)');
    end
end

[~, peakIdx] = max(xc);
peakLagSec = lags(peakIdx) * binSize;

end
