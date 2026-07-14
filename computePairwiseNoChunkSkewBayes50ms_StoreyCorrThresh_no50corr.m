function computePairwiseNoChunkSkewBayes50ms_StoreyCorrThresh_no50corr(alpha, corrThresh, nNullDraws)
% no-chunk pairwise skew analysis comparing H0 against H50

% INPUT:
%   /home/asa7288/pairwise_nochunk_allPairs_ALL_SESSIONS_COMBINED.mat

% OUTPUT:
%   /home/asa7288/pairwiseNoChunkSkewBayes50ms_StoreyCorrThresh.mat

% GOAL:
%   1) Identify significant real interneuron-pyramidal pairs using:
%        - pair-specific circular-shift null p-values
%        - Storey/mafdr q <= alpha
%        - real peak correlation > corrThresh

%   2) Calculate real skew from significant int-pyr peak lags

%   3) Construct an H0 skew distribution by repeatedly sampling the same
%      number of pairs from the full upper-triangle pair pool

%   4) Construct H50 by adding 50 ms to int-pyr peak lags: peakLagH50 = peakLagH0 + 0.050 s for int-pyr pairs only. Int-int and pyr-pyr relative lags remain unchanged.

%   5) Keep peak correlations and correlation nulls unchanged:
%          peakCorrH50 = peakCorrH0
%          nullCorrH50 = nullCorrH0

%      This follows the verified equivalence between physically shifting
%      the interneuron trace by 50 ms and translating the xcorr peak lag
%      by 50 ms.

%   6) Compare the same observed real skew against H0 and H50:
%
%          pH0  = fraction of |H0 skew|  >= |real skew|
%          pH50 = fraction of |H50 skew| >= |real skew|
%          evidence ratio H0/H50 = pH0 / pH50

% NOTES:
%   - no EMG windows
%   - no trial averaging
%   - no firing-rate files are reloaded
%   - no H50 xcorr sweeps are recomputed
%   - no H50 shift-null matrices are recomputed
%   - neuron ordering is [interneurons, pyramidal]

% RUN: computePairwiseNoChunkSkewBayes50ms_StoreyCorrThresh
% or:
%   computePairwiseNoChunkSkewBayes50ms_StoreyCorrThresh(0.05, 0.05, 100)

if nargin < 1 || isempty(alpha)
    alpha = 0.05;
end

if nargin < 2 || isempty(corrThresh)
    corrThresh = 0.05;
end

if nargin < 3 || isempty(nNullDraws)
    nNullDraws = 100;
end

if exist('mafdr', 'file') ~= 2
    error('mafdr was not found. The Bioinformatics Toolbox is required.');
end

%% ---------------- settings ----------------

combinedFile = ...
    '/home/asa7288/pairwise_nochunk_allPairs_ALL_SESSIONS_COMBINED.mat';

outFile = ...
    '/home/asa7288/pairwiseNoChunkSkewBayes50ms_StoreyCorrThresh.mat';

lagShiftSec = 0.050;

rng(0);

%% ---------------- load combined results ----------------

if ~isfile(combinedFile)
    error('Combined input file not found:\n%s', combinedFile);
end

S = load(combinedFile, 'allSessions');
allSessions = S.allSessions;

nSess = numel(allSessions.sessions);

%% ---------------- initialize output ----------------

results = struct();

results.sourceFile = combinedFile;
results.alpha = alpha;
results.corrThresh = corrThresh;
results.nNullDraws = nNullDraws;
results.lagShiftSec = lagShiftSec;

results.analysisNote = [ ...
    'No-chunk pairwise skew analysis. Significant int-pyr pairs were ' ...
    'selected using pair-specific shift-null p-values, Storey/mafdr, ' ...
    'and a minimum peak-correlation threshold. H50 was constructed by ' ...
    'adding 50 ms to int-pyr peak lags while leaving peak correlations ' ...
    'and pair-specific shift-null correlations unchanged.'];

results.sessions = cell(1, nSess);

%% ========================================================================
%% loop through sessions
%% ========================================================================

for sessInd = 1:nSess

    sess = allSessions.sessions(sessInd);

    %% ---------- determine animal ID ----------

    if isfield(sess, 'mouseID') && ~isempty(sess.mouseID)
        animalID = sess.mouseID;
    elseif isfield(sess, 'baseDir') && ~isempty(sess.baseDir)
        animalID = regexp(sess.baseDir, 'D\d+', 'match', 'once');

        if isempty(animalID)
            animalID = sprintf('session_%02d', sessInd);
        end
    else
        animalID = sprintf('session_%02d', sessInd);
    end

    fprintf('\n========================================\n');
    fprintf('processing %s\n', animalID);
    fprintf('========================================\n');

    %% ---------- validate saved session ----------

    requiredFields = { ...
        'peakCorrMatAll', ...
        'peakLagMatAll', ...
        'nullCorrMatAllShifts', ...
        'nInt', ...
        'nPyr', ...
        'nAll'};

    missingRequiredField = false;

    for fieldInd = 1:numel(requiredFields)
        fieldName = requiredFields{fieldInd};

        if ~isfield(sess, fieldName) || isempty(sess.(fieldName))
            warning('%s is missing required field %s.', animalID, fieldName);
            missingRequiredField = true;
        end
    end

    if missingRequiredField
        warning('Skipping %s because required data are missing.', animalID);
        continue;
    end

    if isfield(sess, 'realRowsLoaded') && ...
       ~isempty(sess.realRowsLoaded) && ...
       nnz(sess.realRowsLoaded) < sess.nAll

        warning('%s only has %d/%d real rows loaded.', ...
            animalID, nnz(sess.realRowsLoaded), sess.nAll);
    end

    %% ---------- extract saved H0 data ----------

    nInt = sess.nInt;
    nPyr = sess.nPyr;
    nAll = sess.nAll;

    peakCorrH0 = sess.peakCorrMatAll;
    peakLagH0 = sess.peakLagMatAll;
    nullCorrH0 = sess.nullCorrMatAllShifts;

    if size(peakCorrH0,1) ~= nAll || size(peakCorrH0,2) ~= nAll
        warning('%s peakCorrH0 dimensions do not match nAll. Skipping.', ...
            animalID);
        continue;
    end

    if size(peakLagH0,1) ~= nAll || size(peakLagH0,2) ~= nAll
        warning('%s peakLagH0 dimensions do not match nAll. Skipping.', ...
            animalID);
        continue;
    end

    if size(nullCorrH0,1) ~= nAll || size(nullCorrH0,2) ~= nAll
        warning('%s nullCorrH0 dimensions do not match nAll. Skipping.', ...
            animalID);
        continue;
    end

    fprintf('nInt=%d | nPyr=%d | nAll=%d | shift nulls=%d\n', ...
        nInt, nPyr, nAll, size(nullCorrH0,3));

    %% ====================================================================
    %% real observed int-pyr pair set
    %% ====================================================================

    [actualRows, actualCols] = getIntPyrPairs(nInt, nPyr);

    actual = computePairSetStatsMAFDRSkew_PairSpecificNull( ...
        peakCorrH0, ...
        peakLagH0, ...
        nullCorrH0, ...
        actualRows, ...
        actualCols, ...
        alpha, ...
        corrThresh);

    fprintf(['%s real int-pyr pairs: significant=%d/%d | ' ...
             'valid tests=%d | skew=%.6f\n'], ...
        animalID, ...
        actual.nSigFDR, ...
        actual.nPairsNominal, ...
        actual.nValidTests, ...
        actual.skew);

    if isnan(actual.skew)
        warning(['%s real skew is NaN. This usually means fewer than two ' ...
                 'significant pairs or zero variation among their lags.'], ...
            animalID);
    end

    %% ====================================================================
    %% H0 skew distribution
    %% ====================================================================

    nullSkewH0 = computeNullSkewsFromPairPool( ...
        peakCorrH0, ...
        peakLagH0, ...
        nullCorrH0, ...
        nInt, ...
        nPyr, ...
        nAll, ...
        numel(actualRows), ...
        nNullDraws, ...
        alpha, ...
        corrThresh, ...
        0);

    %% ====================================================================
    %% construct H50 using lag translation only
    %% ====================================================================

    % Correlation magnitude is unchanged by the 50 ms temporal
    % translation. Therefore significance testing uses the same real
    % correlations and the same shift-null correlations.
    peakCorrH50 = peakCorrH0;
    nullCorrH50 = nullCorrH0;

    % Add 50 ms only to the int-pyr pair lags.
    peakLagH50 = addLagToIntPyrPairs( ...
        peakLagH0, ...
        nInt, ...
        nPyr, ...
        lagShiftSec);

    %% ---------- optional range diagnostic ----------

    if isfield(sess, 'lags') && ~isempty(sess.lags) && ...
       isfield(sess, 'binSize') && ~isempty(sess.binSize)

        maxSavedLagSec = max(abs(sess.lags(:))) * sess.binSize;

        intPyrH50Lags = peakLagH50( ...
            sub2ind([nAll nAll], actualRows, actualCols));

        nOutsideOriginalSweep = nnz( ...
            abs(intPyrH50Lags) > maxSavedLagSec);

        if nOutsideOriginalSweep > 0
            warning(['%s has %d H50 int-pyr lags outside the original ' ...
                     '±%.3f s lag-sweep range after adding 50 ms. ' ...
                     'They are retained without clipping.'], ...
                animalID, ...
                nOutsideOriginalSweep, ...
                maxSavedLagSec);
        end
    else
        maxSavedLagSec = NaN;
        nOutsideOriginalSweep = NaN;
    end

    %% ====================================================================
    %% H50 skew distribution
    %% ====================================================================

    nullSkewH50 = computeNullSkewsFromPairPool( ...
        peakCorrH50, ...
        peakLagH50, ...
        nullCorrH50, ...
        nInt, ...
        nPyr, ...
        nAll, ...
        numel(actualRows), ...
        nNullDraws, ...
        alpha, ...
        corrThresh, ...
        lagShiftSec);

    %% ====================================================================
    %% valid skew values
    %% ====================================================================

    validH0 = nullSkewH0( ...
        ~isnan(nullSkewH0) & isfinite(nullSkewH0));

    validH50 = nullSkewH50( ...
        ~isnan(nullSkewH50) & isfinite(nullSkewH50));

    %% ====================================================================
    %% Bayes-style tail probabilities
    %% ====================================================================

    if isnan(actual.skew) || isempty(validH0)
        pValH0 = NaN;
    else
        pValH0 = ...
            (sum(abs(validH0) >= abs(actual.skew)) + 1) / ...
            (numel(validH0) + 1);
    end

    if isnan(actual.skew) || isempty(validH50)
        pValH50 = NaN;
    else
        pValH50 = ...
            (sum(abs(validH50) >= abs(actual.skew)) + 1) / ...
            (numel(validH50) + 1);
    end

    if isnan(pValH0) || isnan(pValH50) || pValH50 == 0
        evidenceRatio_H0_over_H50 = NaN;
    else
        evidenceRatio_H0_over_H50 = pValH0 / pValH50;
    end

    %% ====================================================================
    %% save session output
    %% ====================================================================

    R = struct();

    R.animalID = animalID;

    if isfield(sess, 'baseDir')
        R.baseDir = sess.baseDir;
    else
        R.baseDir = '';
    end

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
    R.nullCorrH0 = nullCorrH0;

    R.peakCorrH50 = peakCorrH50;
    R.peakLagH50 = peakLagH50;
    R.nullCorrH50 = nullCorrH50;

    R.maxSavedLagSec = maxSavedLagSec;
    R.nH50IntPyrLagsOutsideOriginalSweep = nOutsideOriginalSweep;

    R.h50Construction = [ ...
        'peakCorrH50 = peakCorrH0; ' ...
        'nullCorrH50 = nullCorrH0; ' ...
        'peakLagH50 = peakLagH0 + 0.050 s for int-pyr pairs only'];

    results.sessions{sessInd} = R;

    fprintf(['%s skew comparison: real=%.6f | valid H0=%d | ' ...
             'valid H50=%d\n'], ...
        animalID, ...
        actual.skew, ...
        numel(validH0), ...
        numel(validH50));

    fprintf(['%s evidence: pH0=%.6f | pH50=%.6f | ' ...
             'ratio H0/H50=%.6f\n'], ...
        animalID, ...
        pValH0, ...
        pValH50, ...
        evidenceRatio_H0_over_H50);
end

%% ---------------- save all sessions ----------------

save(outFile, 'results', '-v7.3');

fprintf('\nsaved:\n%s\n', outFile);

end

%% ========================================================================
%% helper: real int-pyr pair indices
%% ========================================================================

function [rows, cols] = getIntPyrPairs(nInt, nPyr)
% Ordering:
%   1:nInt              = interneurons
%   nInt+1:nInt+nPyr    = pyramidal neurons

[rr, cc] = ndgrid(1:nInt, nInt + (1:nPyr));

rows = rr(:);
cols = cc(:);

end

%% ========================================================================
%% helper: all upper-triangle pairs
%% ========================================================================

function [rows, cols] = getAllUpperPairs(nAll)

mask = triu(true(nAll), 1);
[rows, cols] = find(mask);

end

%% ========================================================================
%% helper: add 50 ms only to int-pyr lags
%% ========================================================================

function peakLagH50 = addLagToIntPyrPairs( ...
    peakLagH0, nInt, nPyr, lagShiftSec)

nAll = nInt + nPyr;

peakLagH50 = peakLagH0;

[intRows, pyrCols] = getIntPyrPairs(nInt, nPyr);

for pairInd = 1:numel(intRows)

    r = intRows(pairInd);
    c = pyrCols(pairInd);

    if ~isnan(peakLagH0(r,c)) && isfinite(peakLagH0(r,c))
        shiftedLag = peakLagH0(r,c) + lagShiftSec;

        peakLagH50(r,c) = shiftedLag;

        % Preserve the same symmetric-matrix convention used by the
        % combined pairwise files.
        if c <= nAll && r <= nAll
            peakLagH50(c,r) = shiftedLag;
        end
    end
end

end

%% ========================================================================
%% helper: significant-pair statistics with pair-specific nulls
%% ========================================================================

function out = computePairSetStatsMAFDRSkew_PairSpecificNull( ...
    realMat, lagMat, nullMat, rows, cols, alpha, corrThresh, lagSignVec)

if nargin < 8 || isempty(lagSignVec)
    lagSignVec = ones(numel(rows),1);
end

nPairs = numel(rows);

realVals = nan(nPairs,1);
lagVals = nan(nPairs,1);
pVals = nan(nPairs,1);
nNullPerPair = zeros(nPairs,1);

for pairInd = 1:nPairs

    r = rows(pairInd);
    c = cols(pairInd);

    thisReal = realMat(r,c);
    thisLag = lagMat(r,c);

    realVals(pairInd) = thisReal;
    lagVals(pairInd) = lagSignVec(pairInd) * thisLag;

    if ndims(nullMat) ~= 3
        continue;
    end

    thisNull = squeeze(nullMat(r,c,:));

    thisNull = thisNull( ...
        ~isnan(thisNull) & isfinite(thisNull));

    nNullPerPair(pairInd) = numel(thisNull);

    if isnan(thisReal) || ~isfinite(thisReal) || ...
       isnan(thisLag) || ~isfinite(thisLag) || ...
       isempty(thisNull)

        continue;
    end

    % One-sided correlation p-value:
    % how often is shifted-null correlation >= real correlation?
    pVals(pairInd) = ...
        (sum(thisNull >= thisReal) + 1) / ...
        (numel(thisNull) + 1);
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

sigFDR(validP) = ...
    (qVals(validP) <= alpha) & ...
    (realVals(validP) > corrThresh);

sigLagVec = lagVals(sigFDR);

sigLagVec = sigLagVec( ...
    ~isnan(sigLagVec) & isfinite(sigLagVec));

out = struct();

out.rows = rows;
out.cols = cols;

out.realVals = realVals;
out.lagVals = lagVals;
out.pVals = pVals;
out.pFDRVals = pFDRVals;
out.qVals = qVals;
out.nNullPerPair = nNullPerPair;

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

%% ========================================================================
%% helper: matched-size random pair-pool skew draws
%% ========================================================================

function nullSkews = computeNullSkewsFromPairPool( ...
    peakCorr, peakLag, nullCorr, ...
    nInt, nPyr, nAll, nSamplePairs, ...
    nNullDraws, alpha, corrThresh, lagShiftSec)

[poolRows, poolCols] = getAllUpperPairs(nAll);

nPoolPairs = numel(poolRows);

if nSamplePairs > nPoolPairs
    error(['Requested %d sampled pairs, but only %d upper-triangle ' ...
           'pairs are available.'], ...
        nSamplePairs, nPoolPairs);
end

nullSkews = nan(nNullDraws,1);

for drawInd = 1:nNullDraws

    drawIdx = randperm(nPoolPairs, nSamplePairs);

    drawRows = poolRows(drawIdx);
    drawCols = poolCols(drawIdx);

    % Upper-triangle pair orientation is arbitrary for the all-pair null.
    % Randomly reverse lag signs to avoid orientation bias.
    lagSignVec = ones(nSamplePairs,1);

    flipMask = rand(nSamplePairs,1) > 0.5;
    lagSignVec(flipMask) = -1;

    nullDraw = computePairSetStatsMAFDRSkew_PairSpecificNull( ...
        peakCorr, ...
        peakLag, ...
        nullCorr, ...
        drawRows, ...
        drawCols, ...
        alpha, ...
        corrThresh, ...
        lagSignVec);

    nullSkews(drawInd) = nullDraw.skew;
end

end

%% ========================================================================
%% helper: skew metric
%% ========================================================================

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
    skewVal = ...
        (mean(x, 'omitnan') - median(x, 'omitnan')) / sd;
end

end
