function computePairwiseNoChunkSkewBayes50ms_StoreyCorrThresh(alpha, corrThresh, nNullDraws)
% no-chunk pairwise skew analysis comparing H0 against H50

% INPUT DATA: /home/asa7288/pairwise_nochunk_allPairs_ALL_SESSIONS_COMBINED.mat

% GOAL:
%   1) Identify significant real interneuron-pyramidal pairs using:
%        - pair-specific circular-shift null p-values
%        - Storey/mafdr q <= alpha
%        - real peak correlation > corrThresh
%   2) Calculate the real skew from significant int-pyr peak lags
%   3) Construct an H0 skew distribution by repeatedly sampling the same
%      number of pairs from the full upper-triangle pair pool
%   4) Construct H50 pairwise results by shifting interneurons 50 ms
%      earlier relative to pyramidal neurons
%   5) Construct an H50 null-correlation distribution using the same
%      per-neuron circular-shift amounts saved by the Quest shift jobs
%   6) Compare the same real skew against the H0 and H50 skew
%      distributions using:
%         pH0  = fraction of |H0 skew|  >= |real skew|
%         pH50 = fraction of |H50 skew| >= |real skew|
%         evidence ratio H0/H50 = pH0 / pH50

% NOTES:
%   - no EMG windows
%   - no trial averaging
%   - continuous cortex firing-rate traces
%   - neuron ordering is [interneurons, pyramidal]
%   - H50 int-pyr full xcorr curves are derived from the saved H0 curves,
%     avoiding another extremely slow full continuous lag sweep

% RUN: computePairwiseNoChunkSkewBayes50ms_StoreyCorrThresh
% or: computePairwiseNoChunkSkewBayes50ms_StoreyCorrThresh(0.05, 0.05, 100)

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
    'and a minimum real peak-correlation threshold.'];

results.sessions = cell(1, nSess);

%% ========================================================================
%% loop through sessions
%% ========================================================================

for sessInd = 1:nSess

    sess = allSessions.sessions(sessInd);

    if isfield(sess, 'mouseID') && ~isempty(sess.mouseID)
        animalID = sess.mouseID;
    else
        animalID = regexp(sess.baseDir, 'D\d+', 'match', 'once');

        if isempty(animalID)
            animalID = sprintf('session_%02d', sessInd);
        end
    end

    fprintf('\n========================================\n');
    fprintf('processing %s\n', animalID);
    fprintf('========================================\n');

    %% ---------- validate saved session ----------

    if isempty(sess.baseDir) || ...
       isempty(sess.xcMatAll) || ...
       isempty(sess.peakCorrMatAll) || ...
       isempty(sess.peakLagMatAll)

        warning('%s is missing real pairwise results. Skipping.', animalID);
        continue;
    end

    if isempty(sess.nullCorrMatAllShifts)
        warning('%s is missing shifted-null matrices. Skipping.', animalID);
        continue;
    end

    if isempty(sess.shiftAmtPerNeuronAll)
        warning('%s is missing saved per-neuron shift amounts. Skipping.', animalID);
        continue;
    end

    if isfield(sess, 'realRowsLoaded') && ...
       ~isempty(sess.realRowsLoaded) && ...
       nnz(sess.realRowsLoaded) < sess.nAll

        warning('%s only has %d/%d real rows loaded.', ...
            animalID, nnz(sess.realRowsLoaded), sess.nAll);
    end

    nInt = sess.nInt;
    nPyr = sess.nPyr;
    nAll = sess.nAll;

    binSize = sess.binSize;
    lags = sess.lags(:)';

    xcMatH0 = sess.xcMatAll;
    peakCorrH0 = sess.peakCorrMatAll;
    peakLagH0 = sess.peakLagMatAll;

    nullCorrH0 = sess.nullCorrMatAllShifts;

    if size(nullCorrH0,1) ~= nAll || size(nullCorrH0,2) ~= nAll
        warning('%s H0 null matrix dimensions do not match nAll. Skipping.', ...
            animalID);
        continue;
    end

    %% ---------- real int-pyr pair set ----------

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
                 'significant pairs or zero lag variance.'], animalID);
    end

    %% ---------- H0 skew distribution ----------

    nullSkewH0 = computeNullSkewsFromPairPool( ...
        peakCorrH0, ...
        peakLagH0, ...
        nullCorrH0, ...
        nAll, ...
        numel(actualRows), ...
        nNullDraws, ...
        alpha, ...
        corrThresh);

    %% ---------- construct H50 pairwise matrices ----------

    % Rather than rerunning every full continuous pairwise lag sweep,
    % transform the already-saved H0 xcorr curves.
    %
    % For an upper-triangle int-pyr pair:
    %
    %   H50_xc(L) = H0_xc(L + 50 ms)
    %
    % because the first trace is the interneuron trace and it is shifted
    % earlier by 50 ms.

    [xcMatH50, peakCorrH50, peakLagH50] = ...
        deriveH50MatricesFromSavedXCorr( ...
            xcMatH0, ...
            lags, ...
            binSize, ...
            nInt, ...
            nPyr, ...
            lagShiftSec);

    %% ---------- construct H50 continuous traces ----------

    [allFRsH50, okFR] = buildNoChunkH50Traces( ...
        sess.baseDir, ...
        sessInd, ...
        lagShiftSec, ...
        binSize);

    if ~okFR
        warning('Could not construct H50 traces for %s. Skipping.', animalID);
        continue;
    end

    if size(allFRsH50,1) ~= nAll
        warning(['%s H50 trace count=%d, but combined file says nAll=%d. ' ...
                 'Skipping.'], animalID, size(allFRsH50,1), nAll);
        continue;
    end

    %% ---------- construct H50 pair-specific shift null ----------

    % Use the exact same shift amounts saved from the H0 Quest jobs, but
    % apply them to the H50 traces. Each pair therefore still gets one
    % null correlation per shift job.

    nullCorrH50 = computeNoChunkShiftNullsFromSavedShiftAmounts( ...
        allFRsH50, ...
        sess.shiftAmtPerNeuronAll);

    %% ---------- H50 skew distribution ----------

    nullSkewH50 = computeNullSkewsFromPairPool( ...
        peakCorrH50, ...
        peakLagH50, ...
        nullCorrH50, ...
        nAll, ...
        numel(actualRows), ...
        nNullDraws, ...
        alpha, ...
        corrThresh);

    %% ---------- valid skew values ----------

    validH0 = nullSkewH0( ...
        ~isnan(nullSkewH0) & isfinite(nullSkewH0));

    validH50 = nullSkewH50( ...
        ~isnan(nullSkewH50) & isfinite(nullSkewH50));

    %% ---------- Bayes-style tail probabilities ----------

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

    %% ---------- save session output ----------

    R = struct();

    R.animalID = animalID;
    R.baseDir = sess.baseDir;
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

    R.xcMatH0 = xcMatH0;
    R.peakCorrH0 = peakCorrH0;
    R.peakLagH0 = peakLagH0;

    R.xcMatH50 = xcMatH50;
    R.peakCorrH50 = peakCorrH50;
    R.peakLagH50 = peakLagH50;

    R.nullCorrH0 = nullCorrH0;
    R.nullCorrH50 = nullCorrH50;

    results.sessions{sessInd} = R;

    fprintf(['%s skew comparison: real=%.6f | valid H0=%d | ' ...
             'valid H50=%d\n'], ...
        animalID, actual.skew, numel(validH0), numel(validH50));

    fprintf(['%s evidence: pH0=%.6f | pH50=%.6f | ' ...
             'ratio H0/H50=%.6f\n'], ...
        animalID, pValH0, pValH50, evidenceRatio_H0_over_H50);
end

%% ---------------- save all sessions ----------------

save(outFile, 'results', '-v7.3');

fprintf('\nsaved:\n%s\n', outFile);

end

%% ========================================================================
%% helper: real int-pyr pair indices
%% ========================================================================

function [rows, cols] = getIntPyrPairs(nInt, nPyr)
% Ordering in the no-chunk pairwise files is:
%
%   rows 1:nInt          = interneurons
%   rows nInt+1:nAll     = pyramidal neurons
%
% Thus all real int-pyr pairs are in the upper triangle.

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
    % How often is the shifted-null correlation at least as large as real?
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
    peakCorr, peakLag, nullCorr, nAll, nSamplePairs, ...
    nNullDraws, alpha, corrThresh)

[poolRows, poolCols] = getAllUpperPairs(nAll);

nPoolPairs = numel(poolRows);

if nSamplePairs > nPoolPairs
    error(['Requested %d sampled pairs, but only %d upper-triangle ' ...
           'pairs are available.'], nSamplePairs, nPoolPairs);
end

nullSkews = nan(nNullDraws,1);

for drawInd = 1:nNullDraws

    drawIdx = randperm(nPoolPairs, nSamplePairs);

    drawRows = poolRows(drawIdx);
    drawCols = poolCols(drawIdx);

    % Since upper-triangle pairs have an arbitrary first/second
    % orientation, randomly reverse the lag signs to avoid introducing an
    % orientation bias into the sampled all-pair null.
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
%% helper: derive H50 xcorr and peak matrices from saved H0 xcorr
%% ========================================================================

function [xcMatH50, peakCorrH50, peakLagH50] = ...
    deriveH50MatricesFromSavedXCorr( ...
        xcMatH0, lags, binSize, nInt, nPyr, lagShiftSec)

nAll = nInt + nPyr;
nLags = numel(lags);

lagShiftBins = round(lagShiftSec / binSize);

xcMatH50 = nan(size(xcMatH0));
peakCorrH50 = nan(nAll,nAll);
peakLagH50 = nan(nAll,nAll);

for i = 1:nAll
    for j = (i+1):nAll

        xcH0 = squeeze(xcMatH0(i,j,:))';

        if isempty(xcH0) || all(isnan(xcH0))
            continue;
        end

        isIntPyrPair = (i <= nInt) && (j > nInt);

        if isIntPyrPair
            % H50_xc(L) = H0_xc(L + 50 ms).
            %
            % Therefore target lag index t receives source index
            % t + lagShiftBins.
            xcH50 = nan(1,nLags);

            validTarget = 1:(nLags - lagShiftBins);
            validSource = (1 + lagShiftBins):nLags;

            xcH50(validTarget) = xcH0(validSource);
        else
            % Shifting both interneurons equally leaves int-int relative
            % timing unchanged. Pyr-pyr traces are not shifted.
            xcH50 = xcH0;
        end

        xcMatH50(i,j,:) = xcH50;
        xcMatH50(j,i,:) = xcH50;

        validXC = ~isnan(xcH50) & isfinite(xcH50);

        if any(validXC)
            validIndices = find(validXC);
            [peakCorr, localPeakIndex] = max(xcH50(validXC));
            peakIndex = validIndices(localPeakIndex);

            peakLag = lags(peakIndex) * binSize;

            peakCorrH50(i,j) = peakCorr;
            peakCorrH50(j,i) = peakCorr;

            peakLagH50(i,j) = peakLag;
            peakLagH50(j,i) = peakLag;
        end
    end
end

end

%% ========================================================================
%% helper: load continuous traces and shift interneurons earlier by 50 ms
%% ========================================================================

function [allFRsH50, ok] = buildNoChunkH50Traces( ...
    baseDir, sessInd, lagShiftSec, binSize)

ok = false;
allFRsH50 = [];

frFile = fullfile( ...
    baseDir, ...
    'NeuralFiringRates1msBins10msGauss.mat');

clsFile = ...
    '/home/asa7288/Transfer/AA_classifications.mat';

if ~isfile(frFile)
    warning('Missing firing-rate file:\n%s', frFile);
    return;
end

if ~isfile(clsFile)
    warning('Missing classification file:\n%s', clsFile);
    return;
end

F = load(frFile, 'cortexFRs', 'cortexInds');
C = load(clsFile, 'classifications');

regionClass = C.classifications{sessInd,1}(F.cortexInds);

interFRs = F.cortexFRs(regionClass == 1,:);
pyrFRs = F.cortexFRs(regionClass == 0,:);

lagShiftBins = round(lagShiftSec / binSize);

if lagShiftBins < 1 || lagShiftBins >= size(interFRs,2)
    warning('Invalid H50 shift of %d bins.', lagShiftBins);
    return;
end

interFRsH50 = nan(size(interFRs));

% Shift interneuron activity earlier by 50 ms:
%
% new(t) = original(t + 50 ms)
interFRsH50(:,1:end-lagShiftBins) = ...
    interFRs(:,1+lagShiftBins:end);

allFRsH50 = [interFRsH50; pyrFRs];

ok = true;

end

%% ========================================================================
%% helper: H50 pair-specific shift-null matrices
%% ========================================================================

function nullCorrMatAllShifts = ...
    computeNoChunkShiftNullsFromSavedShiftAmounts( ...
        allFRs, shiftAmtPerNeuronAll)

nAll = size(allFRs,1);
nShifts = size(shiftAmtPerNeuronAll,2);

if size(shiftAmtPerNeuronAll,1) ~= nAll
    error(['shiftAmtPerNeuronAll has %d rows, but allFRs has %d ' ...
           'neurons.'], size(shiftAmtPerNeuronAll,1), nAll);
end

nullCorrMatAllShifts = nan(nAll,nAll,nShifts);

for shiftInd = 1:nShifts

    shiftAmtPerNeuron = shiftAmtPerNeuronAll(:,shiftInd);

    if any(isnan(shiftAmtPerNeuron))
        warning('Shift %d has missing per-neuron shift values.', shiftInd);
        continue;
    end

    allFRsShifted = allFRs;

    for neuronInd = 1:nAll
        allFRsShifted(neuronInd,:) = circshift( ...
            allFRs(neuronInd,:), ...
            [0 shiftAmtPerNeuron(neuronInd)]);
    end

    nullCorrThisShift = nan(nAll,nAll);

    for i = 1:nAll

        unshiftedTrace = allFRs(i,:);

        for j = (i+1):nAll

            shiftedTrace = allFRsShifted(j,:);

            valid = ...
                ~isnan(unshiftedTrace) & ...
                ~isnan(shiftedTrace);

            if nnz(valid) > 2
                r = corr( ...
                    unshiftedTrace(valid)', ...
                    shiftedTrace(valid)');

                nullCorrThisShift(i,j) = r;
                nullCorrThisShift(j,i) = r;
            end
        end
    end

    nullCorrMatAllShifts(:,:,shiftInd) = nullCorrThisShift;

    if mod(shiftInd,10) == 0 || shiftInd == nShifts
        fprintf('  rebuilt H50 shift null %d/%d\n', ...
            shiftInd, nShifts);
    end
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
