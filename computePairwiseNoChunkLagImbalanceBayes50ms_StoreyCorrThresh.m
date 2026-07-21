function computePairwiseNoChunkLagImbalanceBayes50ms_StoreyCorrThresh(alpha, corrThresh, nNullDraws, maxPermutationTries)
% no-chunk pairwise lag-imbalance analysis comparing:

%   H0 = no additional lag shift
%   H+50 = add +50 ms to every selected, oriented lag
%   H-50 = subtract 50 ms from every selected, oriented lag

% INPUT: /home/asa7288/pairwise_nochunk_allPairs_ALL_SESSIONS_COMBINED.mat
% OUTPUT: /home/asa7288/pairwiseNoChunkLagImbalanceBayes50ms_StoreyCorrThresh.mat

% CONFIRMED PERMUTATION LOGIC:
%   1) Determine significant upper-triangle pairs once using:
%        - pair-specific circular-shift correlation nulls
%        - Storey/mafdr q <= alpha
%        - real peak correlation > corrThresh
%   2) Calculate the actual lag imbalance using significant true
%      interneuron-pyramidal pairs, oriented as real int -> real pyr:
%         lagImbalance = (nPositive - nNegative) / (nPositive + nNegative)
%      Exact zero lags are excluded from the denominator.
%   3) For every permutation draw:
%        - randomly choose nInt neurons as pseudo-interneurons
%        - assign the remaining nPyr neurons as pseudo-pyramidal
%        - form all pseudo-int x pseudo-pyr pairs
%        - retain only pairs that passed the fixed upper-triangle
%          Storey + correlation-threshold significance test
%        - orient each lag as pseudo-int -> pseudo-pyr:
%              if pseudoInt < pseudoPyr
%                  lag = peakLagMat(pseudoInt, pseudoPyr)
%              else
%                  lag = -peakLagMat(pseudoPyr, pseudoInt)
%              end

%        - if fewer significant eligible pairs are available than num significant real int-pyr pairs, reject the label permutation and redraw it
%        - otherwise randomly select exactly the same number of pairs as the number of significant real int-pyr pairs
%        - apply the model shift:
%              H0: +0.000 s
%              H+50: +0.050 s
%              H-50: -0.050 s
%        - calculate one lag-imbalance value

%   4) H0, H+50, and H-50 are generated using separate, independent sets of neuron-label permutations.

%   5) Compare the fixed actual lag imbalance against each null distribution using empirical two-sided tail probabilities:
%         evidence ratio H0/H+50 = pH0 / pH50
%         evidence ratio H0/H-50 = pH0 / pHneg50

% RUN: computePairwiseNoChunkLagImbalanceBayes50ms_StoreyCorrThresh
% or: computePairwiseNoChunkLagImbalanceBayes50ms_StoreyCorrThresh(0.05, 0.05, 100, 1000)

if nargin < 1 || isempty(alpha)
    alpha = 0.05;
end

if nargin < 2 || isempty(corrThresh)
    corrThresh = 0.05;
end

if nargin < 3 || isempty(nNullDraws)
    nNullDraws = 100;
end

if nargin < 4 || isempty(maxPermutationTries)
    maxPermutationTries = 1000;
end

if exist('mafdr', 'file') ~= 2
    error(['mafdr was not found. The Bioinformatics Toolbox must be ' 'available on the MATLAB path.']);
end

%% ---------------- settings ----------------

combinedFile = '/home/asa7288/pairwise_nochunk_allPairs_ALL_SESSIONS_COMBINED.mat';

outFile = '/home/asa7288/pairwiseNoChunkLagImbalanceBayes50ms_StoreyCorrThresh.mat';

lagShiftSec = 0.050;

% Reproducible overall run.
% H0, H+50, and H-50 still use independent draws because each helper call
% continues advancing the random-number stream.
rng(0);

%% ---------------- load combined results ----------------

if ~isfile(combinedFile)
    error('Combined input file not found:\n%s', combinedFile);
end

S = load(combinedFile, 'allSessions');

if ~isfield(S, 'allSessions') || ~isfield(S.allSessions, 'sessions')

    error('Input file does not contain allSessions.sessions.');
end

allSessions = S.allSessions;
nSess = numel(allSessions.sessions);

%% ---------------- initialize output ----------------

results = struct();
results.sourceFile = combinedFile;
results.outputFile = outFile;
results.alpha = alpha;
results.corrThresh = corrThresh;
results.nNullDraws = nNullDraws;
results.maxPermutationTries = maxPermutationTries;
results.lagShiftSec = lagShiftSec;
results.fdrMethod = 'Storey mafdr across all valid upper-triangle pair p-values';
results.significanceRule = 'q <= alpha and real peak correlation > corrThresh';
results.lagImbalanceDefinition = '(nPositive - nNegative) / (nPositive + nNegative); exact zero lags excluded';
results.permutationMethod = ['Randomly reassign neuron-type labels while preserving nInt and nPyr; ' ...
    'form pseudo-int x pseudo-pyr pairs; retain fixed globally significant ' ...
    'upper-triangle pairs; orient lags as pseudo-int to pseudo-pyr; redraw ' ...
    'label assignment when too few eligible pairs are available.'];

results.modelDefinitions = struct('H0_shiftSec', 0, 'H50_shiftSec', +lagShiftSec, 'Hneg50_shiftSec', -lagShiftSec);

results.sessions = cell(1, nSess);

%% ========================================================================
%% loop through sessions
%% ========================================================================

for sessInd = 1:nSess

    sess = allSessions.sessions(sessInd);

    %% ---------- determine animal ID ----------

    if isfield(sess, 'mouseID') && ~isempty(sess.mouseID)
        animalID = char(string(sess.mouseID));

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

    %% ---------- validate required session fields ----------

    requiredFields = {'peakCorrMatAll', 'peakLagMatAll', 'nullCorrMatAllShifts', 'nInt', 'nPyr', 'nAll'};

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

    %% ---------- extract data ----------

    nInt = double(sess.nInt);
    nPyr = double(sess.nPyr);
    nAll = double(sess.nAll);

    peakCorrMat = sess.peakCorrMatAll;
    peakLagMat = sess.peakLagMatAll;
    nullCorrMat = sess.nullCorrMatAllShifts;

    if nInt + nPyr ~= nAll
        warning(['%s has nInt + nPyr = %d, but nAll = %d. ' 'Skipping session.'], animalID, nInt + nPyr, nAll);
        continue;
    end

    if size(peakCorrMat,1) ~= nAll || size(peakCorrMat,2) ~= nAll
        warning('%s peakCorrMatAll dimensions do not match nAll.', animalID);
        continue;
    end

    if size(peakLagMat,1) ~= nAll || size(peakLagMat,2) ~= nAll
        warning('%s peakLagMatAll dimensions do not match nAll.', animalID);
        continue;
    end

    if size(nullCorrMat,1) ~= nAll || size(nullCorrMat,2) ~= nAll

        warning(['%s nullCorrMatAllShifts dimensions do not match ' 'nAll.'], animalID);
        continue;
    end

    fprintf('nInt=%d | nPyr=%d | nAll=%d | shift nulls=%d\n', nInt, nPyr, nAll, size(nullCorrMat,3));

    if isfield(sess, 'realRowsLoaded') && ~isempty(sess.realRowsLoaded) && nnz(sess.realRowsLoaded) < nAll

        warning('%s only has %d/%d real rows loaded.', animalID, nnz(sess.realRowsLoaded), nAll);
    end

    %% ====================================================================
    %% determine significance once across full upper triangle
    %% ====================================================================

    [upperRows, upperCols] = getAllUpperPairs(nAll);

    allPairStats = computePairSignificance_PairSpecificNull(peakCorrMat, peakLagMat, nullCorrMat, upperRows, upperCols, alpha, corrThresh);

    % Fixed matrix indicating which saved upper-triangle pairs passed Storey + correlation threshold
    sigUpperMask = false(nAll, nAll);

    sigLinearInds = sub2ind([nAll nAll], upperRows(allPairStats.sigFDRMask), upperCols(allPairStats.sigFDRMask));
    sigUpperMask(sigLinearInds) = true;
    nSignificantUpperPairs = nnz(allPairStats.sigFDRMask);

    fprintf(['all upper-triangle pairs: significant=%d/%d | ' 'valid tests=%d\n'], nSignificantUpperPairs, allPairStats.nPairsNominal, allPairStats.nValidTests);

    if nSignificantUpperPairs == 0
        warning('%s has no significant upper-triangle pairs. Skipping.', animalID);
        continue;
    end

    %% ====================================================================
    %% calculate actual real int-pyr lag imbalance
    %% ====================================================================

    % True neuron ordering in the saved matrices:
    %   1:nInt = real interneurons
    %   nInt+1:nAll = real pyramidal neurons
    
    % Therefore every real int-pyr pair naturally appears as row < column, and its stored lag already has the desired int -> pyr orientation.

    [actualRows, actualCols] = getIntPyrPairs(nInt, nPyr);

    actualUpperInds = sub2ind([nAll nAll], actualRows, actualCols);

    actualSigMask = sigUpperMask(actualUpperInds);
    actualSigRows = actualRows(actualSigMask);
    actualSigCols = actualCols(actualSigMask);
    actualSigLags = peakLagMat(sub2ind([nAll nAll], actualSigRows, actualSigCols));
    actualSigLags = actualSigLags(:);
    actualSigLags = actualSigLags(isfinite(actualSigLags));

    nActualSigIntPyr = numel(actualSigLags);

    if nActualSigIntPyr == 0
        warning('%s has no significant real int-pyr pairs. Skipping.', animalID);
        continue;
    end

    actualLagCounts = countLagSigns(actualSigLags);

    actualLagImbalance = computeLagImbalance(actualSigLags);

    fprintf(['real int-pyr significant pairs=%d/%d | ' 'negative=%d | zero=%d | positive=%d | ' 'lag imbalance=%.6f\n'], nActualSigIntPyr, ...
        numel(actualRows), actualLagCounts.nNegative, actualLagCounts.nZero, actualLagCounts.nPositive, actualLagImbalance);

    if ~isfinite(actualLagImbalance)
        warning(['%s actual lag imbalance is undefined because all ' 'significant real int-pyr lags were exactly zero.'], animalID);
        continue;
    end

    %% ====================================================================
    %% construct independent H0, H+50, and H-50 distributions
    %% ====================================================================

    fprintf('\nbuilding H0 permutation distribution...\n');

    [nullLagImbalanceH0, permutationTriesH0, selectedLagCountsH0, selectedPositiveCountsH0, selectedNegativeCountsH0, selectedZeroCountsH0, matchedTargetH0, bestEligibleCountsH0] = ...
        computeLabelPermutationLagImbalanceDistribution(peakLagMat, sigUpperMask, nInt, nPyr, nActualSigIntPyr, nNullDraws, 0, maxPermutationTries);

    fprintf('\nbuilding H+50 permutation distribution...\n');
    [nullLagImbalanceH50, permutationTriesH50, selectedLagCountsH50, selectedPositiveCountsH50, selectedNegativeCountsH50, selectedZeroCountsH50, matchedTargetH50, bestEligibleCountsH50] = ...
        computeLabelPermutationLagImbalanceDistribution(peakLagMat, sigUpperMask, nInt, nPyr, nActualSigIntPyr, nNullDraws, +lagShiftSec, maxPermutationTries);

    fprintf('\nbuilding H-50 permutation distribution...\n');

    [nullLagImbalanceHneg50, permutationTriesHneg50, selectedLagCountsHneg50, selectedPositiveCountsHneg50, selectedNegativeCountsHneg50, ...
     selectedZeroCountsHneg50, matchedTargetHneg50, bestEligibleCountsHneg50] = computeLabelPermutationLagImbalanceDistribution(peakLagMat, sigUpperMask, nInt, nPyr, nActualSigIntPyr, ...
            nNullDraws, -lagShiftSec, maxPermutationTries);

    %% ====================================================================
    %% retain valid null values
    %% ====================================================================

    validH0 = nullLagImbalanceH0(isfinite(nullLagImbalanceH0));
    validH50 = nullLagImbalanceH50(isfinite(nullLagImbalanceH50));
    validHneg50 = nullLagImbalanceHneg50(isfinite(nullLagImbalanceHneg50));

    %% ====================================================================
    %% empirical two-sided tail probabilities
    %% ====================================================================

    pValH0 = empiricalTwoSidedTailProbability(actualLagImbalance, validH0);
    pValH50 = empiricalTwoSidedTailProbability(actualLagImbalance, validH50);
    pValHneg50 = empiricalTwoSidedTailProbability(actualLagImbalance, validHneg50);

    %% ====================================================================
    %% Bayes-style evidence ratios
    %% ====================================================================

    evidenceRatio_H0_over_H50 = safeRatio(pValH0, pValH50);

    evidenceRatio_H0_over_Hneg50 = safeRatio(pValH0, pValHneg50);

    %% ====================================================================
    %% calculate descriptive null intervals
    %% ====================================================================

    if isempty(validH0)
        nullCIH0 = [NaN NaN];
    else
        nullCIH0 = prctile(validH0, [2.5 97.5]);
    end

    if isempty(validH50)
        nullCIH50 = [NaN NaN];
    else
        nullCIH50 = prctile(validH50, [2.5 97.5]);
    end

    if isempty(validHneg50)
        nullCIHneg50 = [NaN NaN];
    else
        nullCIHneg50 = prctile(validHneg50, [2.5 97.5]);
    end

    %% ====================================================================
    %% save session output
    %% ====================================================================

    R = struct();

    R.animalID = animalID;
    R.sessInd = sessInd;

    if isfield(sess, 'baseDir')
        R.baseDir = sess.baseDir;
    else
        R.baseDir = '';
    end

    R.nInt = nInt;
    R.nPyr = nPyr;
    R.nAll = nAll;

    R.alpha = alpha;
    R.corrThresh = corrThresh;
    R.nNullDraws = nNullDraws;
    R.maxPermutationTries = maxPermutationTries;
    R.lagShiftSec = lagShiftSec;

    % Fixed upper-triangle significance results
    R.allPairStats = allPairStats;
    R.upperRows = upperRows;
    R.upperCols = upperCols;
    R.sigUpperMask = sigUpperMask;
    R.nSignificantUpperPairs = nSignificantUpperPairs;

    % Actual true int-pyr result
    R.actualRows = actualRows;
    R.actualCols = actualCols;
    R.actualSigMask = actualSigMask;
    R.actualSigRows = actualSigRows;
    R.actualSigCols = actualSigCols;
    R.actualSigLags = actualSigLags;
    R.nActualSigIntPyr = nActualSigIntPyr;
    R.actualLagCounts = actualLagCounts;
    R.actualLagImbalance = actualLagImbalance;

    % Null distributions
    R.nullLagImbalanceH0 = nullLagImbalanceH0;
    R.nullLagImbalanceH50 = nullLagImbalanceH50;
    R.nullLagImbalanceHneg50 = nullLagImbalanceHneg50;

    R.validNullLagImbalanceH0 = validH0;
    R.validNullLagImbalanceH50 = validH50;
    R.validNullLagImbalanceHneg50 = validHneg50;

    R.nullCIH0 = nullCIH0;
    R.nullCIH50 = nullCIH50;
    R.nullCIHneg50 = nullCIHneg50;

    % Number of label-permutation attempts required for each accepted draw
    R.permutationTriesH0 = permutationTriesH0;
    R.permutationTriesH50 = permutationTriesH50;
    R.permutationTriesHneg50 = permutationTriesHneg50;

    % These equal nActualSigIntPyr when matchedTarget is true; fallback draws may use fewer pairs
    R.selectedLagCountsH0 = selectedLagCountsH0;
    R.selectedLagCountsH50 = selectedLagCountsH50;
    R.selectedLagCountsHneg50 = selectedLagCountsHneg50;

    % Whether each draw matched the real significant-pair count exactly
    R.matchedTargetH0 = matchedTargetH0;
    R.matchedTargetH50 = matchedTargetH50;
    R.matchedTargetHneg50 = matchedTargetHneg50;

    % Largest eligible significant-pair count found within the search
    R.bestEligibleCountsH0 = bestEligibleCountsH0;
    R.bestEligibleCountsH50 = bestEligibleCountsH50;
    R.bestEligibleCountsHneg50 = bestEligibleCountsHneg50;

    % Sign counts after applying each model shift
    R.selectedPositiveCountsH0 = selectedPositiveCountsH0;
    R.selectedNegativeCountsH0 = selectedNegativeCountsH0;
    R.selectedZeroCountsH0 = selectedZeroCountsH0;

    R.selectedPositiveCountsH50 = selectedPositiveCountsH50;
    R.selectedNegativeCountsH50 = selectedNegativeCountsH50;
    R.selectedZeroCountsH50 = selectedZeroCountsH50;

    R.selectedPositiveCountsHneg50 = selectedPositiveCountsHneg50;
    R.selectedNegativeCountsHneg50 = selectedNegativeCountsHneg50;
    R.selectedZeroCountsHneg50 = selectedZeroCountsHneg50;

    % Tail probabilities and evidence ratios
    R.pValH0 = pValH0;
    R.pValH50 = pValH50;
    R.pValHneg50 = pValHneg50;

    R.evidenceRatio_H0_over_H50 = evidenceRatio_H0_over_H50;
    R.evidenceRatio_H0_over_Hneg50 = evidenceRatio_H0_over_Hneg50;
    R.lagImbalanceDefinition = '(nPositive - nNegative) / (nPositive + nNegative)';

    R.modelConstruction = [ ...
        'Significance is fixed once across all upper-triangle pairs. ' ...
        'Each null draw randomly reassigns nInt pseudo-interneuron labels ' ...
        'and nPyr pseudo-pyramidal labels. Pseudo-int x pseudo-pyr lags ' ...
        'are oriented deterministically from matrix ordering. ' ...
        'For each draw, up to maxPermutationTries label assignments are tested. ' ...
        'If a permutation reaches nActualSigIntPyr eligible significant pairs, ' ...
        'exactly that many are sampled without replacement. Otherwise, all ' ...
        'eligible pairs from the permutation with the largest eligible count ' ...
        'are used as a fallback. H0 applies ' ...
        '0 ms, H+50 applies +50 ms, and H-50 applies -50 ms after lag ' ...
        'orientation.'];

    results.sessions{sessInd} = R;

    %% ---------- print session summary ----------

    fprintf('\n%s lag-imbalance summary:\n', animalID);
    fprintf('  actual lag imbalance: %.6f\n', actualLagImbalance);
    fprintf('  valid null draws: H0=%d | H+50=%d | H-50=%d\n', numel(validH0), numel(validH50), numel(validHneg50));
    fprintf('  p-values: H0=%.6f | H+50=%.6f | H-50=%.6f\n', pValH0, pValH50, pValHneg50);
    fprintf('  evidence ratios: H0/H+50=%.6f | H0/H-50=%.6f\n', evidenceRatio_H0_over_H50, evidenceRatio_H0_over_Hneg50);
    fprintf(['  median permutation tries: H0=%.1f | H+50=%.1f | ' 'H-50=%.1f\n'], median(permutationTriesH0), ...
        median(permutationTriesH50), median(permutationTriesHneg50));
    fprintf('  exact pair-count matches: H0=%d/%d | H+50=%d/%d | H-50=%d/%d\n', ...
        nnz(matchedTargetH0), nNullDraws, nnz(matchedTargetH50), nNullDraws, ...
        nnz(matchedTargetHneg50), nNullDraws);
end

%% ---------------- save all sessions ----------------

save(outFile, 'results', '-v7.3');

fprintf('\n========================================\n');
fprintf('saved:\n%s\n', outFile);
fprintf('========================================\n');

end

%% ========================================================================
%% helper: true int-pyr pair indices
%% ========================================================================

function [rows, cols] = getIntPyrPairs(nInt, nPyr)
% Saved neuron ordering:
%   1:nInt = true interneurons
%   nInt+1:nInt+nPyr = true pyramidal neurons
%
% Therefore all true int-pyr pairs are stored in the upper triangle and
% already have the desired int -> pyr orientation.

[rr, cc] = ndgrid(1:nInt, nInt + (1:nPyr));

rows = rr(:);
cols = cc(:);

end

%% ========================================================================
%% helper: all saved upper-triangle pair indices
%% ========================================================================

function [rows, cols] = getAllUpperPairs(nAll)

upperMask = triu(true(nAll), 1);
[rows, cols] = find(upperMask);

end

%% ========================================================================
%% helper: calculate fixed Storey + correlation-threshold significance
%% ========================================================================

function out = computePairSignificance_PairSpecificNull(realMat, lagMat, nullMat, rows, cols, alpha, corrThresh)

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
    lagVals(pairInd) = thisLag;

    if ndims(nullMat) ~= 3
        continue;
    end

    thisNull = squeeze(nullMat(r,c,:));
    thisNull = thisNull(isfinite(thisNull));

    nNullPerPair(pairInd) = numel(thisNull);

    if ~isfinite(thisReal) || ~isfinite(thisLag) || isempty(thisNull)

        continue;
    end

    % One-sided empirical correlation p-value:
    % how often is shifted-null correlation at least as large as the
    % observed real peak correlation?
    pVals(pairInd) = (sum(thisNull >= thisReal) + 1) / (numel(thisNull) + 1);
end

validP = isfinite(pVals);

sigUncorrMask = false(size(pVals));
sigUncorrMask(validP) = pVals(validP) <= alpha;

qVals = nan(size(pVals));
pFDRVals = nan(size(pVals));

if any(validP)
    [pFDRtmp, qTmp] = mafdr(pVals(validP));
    pFDRVals(validP) = pFDRtmp;
    qVals(validP) = qTmp;
end

sigFDRMask = false(size(pVals));

sigFDRMask(validP) = (qVals(validP) <= alpha) & (realVals(validP) > corrThresh);

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
out.nSigUncorr = nnz(sigUncorrMask);
out.nSigFDR = nnz(sigFDRMask);
out.sigUncorrMask = sigUncorrMask;
out.sigFDRMask = sigFDRMask;

end

%% ========================================================================
%% helper: neuron-label permutation lag-imbalance distribution
%% ========================================================================

function [nullLagImbalance, permutationTries, selectedLagCounts, selectedPositiveCounts, selectedNegativeCounts, ...
    selectedZeroCounts, matchedTarget, bestEligibleCounts] = computeLabelPermutationLagImbalanceDistribution(peakLagMat, sigUpperMask, nInt, nPyr, ...
        nTargetPairs, nNullDraws, modelShiftSec, maxPermutationTries)

% For every null draw:
%   1) Randomly assign nInt neurons as pseudo-interneurons.
%   2) Assign the remaining nPyr neurons as pseudo-pyramidal.
%   3) Form all pseudo-int x pseudo-pyr pairs.
%   4) Retain pairs that passed the fixed upper-triangle significance mask.
%   5) Orient each lag as pseudo-int -> pseudo-pyr.
%   6) Try up to maxPermutationTries to find a label permutation with at
%      least nTargetPairs eligible significant pairs.
%   7) If one is found, randomly select exactly nTargetPairs without
%      replacement.
%   8) If none is found, use all eligible lags from the permutation with
%      the largest eligible significant-pair count.
%   9) Add modelShiftSec and calculate lag imbalance.

nAll = nInt + nPyr;

nullLagImbalance = nan(nNullDraws,1);
permutationTries = nan(nNullDraws,1);
selectedLagCounts = nan(nNullDraws,1);
selectedPositiveCounts = nan(nNullDraws,1);
selectedNegativeCounts = nan(nNullDraws,1);
selectedZeroCounts = nan(nNullDraws,1);
matchedTarget = false(nNullDraws,1);
bestEligibleCounts = nan(nNullDraws,1);

for drawInd = 1:nNullDraws

    acceptedExactMatch = false;
    bestEligibleCount = -1;
    bestEligibleLags = [];
    bestTry = NaN;
    selectedLags = [];

    for tryCount = 1:maxPermutationTries

        %% ---------- randomly reassign neuron-type labels ----------

        permutedNeuronOrder = randperm(nAll);
        pseudoIntInds = permutedNeuronOrder(1:nInt);
        pseudoPyrInds = permutedNeuronOrder(nInt+1:end);

        %% ---------- gather significant pseudo-int x pseudo-pyr lags ----------

        maxPossiblePairs = nInt * nPyr;
        eligibleLags = nan(maxPossiblePairs,1);
        nEligible = 0;

        for intInd = 1:nInt

            pseudoIntNeuron = pseudoIntInds(intInd);

            for pyrInd = 1:nPyr

                pseudoPyrNeuron = pseudoPyrInds(pyrInd);

                if pseudoIntNeuron < pseudoPyrNeuron
                    upperRow = pseudoIntNeuron;
                    upperCol = pseudoPyrNeuron;
                    orientationSign = 1;
                else
                    upperRow = pseudoPyrNeuron;
                    upperCol = pseudoIntNeuron;
                    orientationSign = -1;
                end

                if ~sigUpperMask(upperRow, upperCol)
                    continue;
                end

                storedLag = peakLagMat(upperRow, upperCol);

                if ~isfinite(storedLag)
                    continue;
                end

                nEligible = nEligible + 1;
                eligibleLags(nEligible) = orientationSign * storedLag;
            end
        end

        eligibleLags = eligibleLags(1:nEligible);

        %% ---------- remember the best available permutation ----------

        if nEligible > bestEligibleCount
            bestEligibleCount = nEligible;
            bestEligibleLags = eligibleLags;
            bestTry = tryCount;
        end

        %% ---------- preferred case: exactly match real pair count ----------

        if nEligible >= nTargetPairs
            chosenIndices = randperm(nEligible, nTargetPairs);
            selectedLags = eligibleLags(chosenIndices);
            acceptedExactMatch = true;
            permutationTries(drawInd) = tryCount;
            matchedTarget(drawInd) = true;
            break;
        end
    end

    %% ---------- fallback: use best available permutation ----------

    if ~acceptedExactMatch

        if isempty(bestEligibleLags)
            error(['No eligible pseudo-int x pseudo-pyr pairs were found ' ...
                   'for draw %d after %d label permutations.'], ...
                   drawInd, maxPermutationTries);
        end

        selectedLags = bestEligibleLags;
        permutationTries(drawInd) = maxPermutationTries;
        matchedTarget(drawInd) = false;

        fprintf(['  fallback for draw %d: target=%d, best available=%d ' ...
                 '(found on try %d)\n'], ...
                 drawInd, nTargetPairs, bestEligibleCount, bestTry);
    end

    bestEligibleCounts(drawInd) = bestEligibleCount;

    %% ---------- apply model shift after lag orientation ----------

    modelLags = selectedLags + modelShiftSec;

    %% ---------- calculate lag imbalance ----------

    thisLagImbalance = computeLagImbalance(modelLags);

    if ~isfinite(thisLagImbalance)
        error(['Lag imbalance was undefined for draw %d because all ' ...
               'selected model lags were exactly zero.'], drawInd);
    end

    nullLagImbalance(drawInd) = thisLagImbalance;

    thisCounts = countLagSigns(modelLags);

    selectedLagCounts(drawInd) = numel(modelLags);
    selectedPositiveCounts(drawInd) = thisCounts.nPositive;
    selectedNegativeCounts(drawInd) = thisCounts.nNegative;
    selectedZeroCounts(drawInd) = thisCounts.nZero;

    if mod(drawInd,10) == 0 || drawInd == nNullDraws
        fprintf(['  completed draw %d/%d | model shift=%+.3f s | ' ...
                 'selected pairs=%d | exact match=%d\n'], ...
                 drawInd, nNullDraws, modelShiftSec, ...
                 selectedLagCounts(drawInd), matchedTarget(drawInd));
    end
end

end

%% ========================================================================
%% helper: lag imbalance
%% ========================================================================

function lagImbalance = computeLagImbalance(lags)
% Exact zero lags are excluded from the denominator.

% Range:
%   -1 = all nonzero lags are negative
%    0 = equal numbers of positive and negative lags
%   +1 = all nonzero lags are positive

lags = lags(isfinite(lags));

nPositive = nnz(lags > 0);
nNegative = nnz(lags < 0);

denominator = nPositive + nNegative;

if denominator == 0
    lagImbalance = NaN;
else
    lagImbalance = ...
        (nPositive - nNegative) / denominator;
end

end

%% ========================================================================
%% helper: count positive, negative, and zero lags
%% ========================================================================

function counts = countLagSigns(lags)

lags = lags(isfinite(lags));

counts = struct();

counts.nNegative = nnz(lags < 0);
counts.nZero = nnz(lags == 0);
counts.nPositive = nnz(lags > 0);

counts.nNonzero = counts.nNegative + counts.nPositive;

end

%% ========================================================================
%% helper: empirical two-sided tail probability
%% ========================================================================

function pVal = empiricalTwoSidedTailProbability(actualValue, nullValues)

nullValues = nullValues(isfinite(nullValues));

if ~isfinite(actualValue) || isempty(nullValues)
    pVal = NaN;
    return;
end

pVal = (sum(abs(nullValues) >= abs(actualValue)) + 1) / (numel(nullValues) + 1);

end

%% ========================================================================
%% helper: safe evidence-ratio calculation
%% ========================================================================

function ratio = safeRatio(numerator, denominator)

if ~isfinite(numerator) || ...
   ~isfinite(denominator) || ...
   denominator == 0

    ratio = NaN;

else

    ratio = numerator / denominator;
end

end
