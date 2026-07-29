function computePairwiseNoChunkLagImbalanceBayes50ms_StoreyCorrThresh_final(alpha, corrThresh, nNullDraws)
% No-chunk pairwise lag-imbalance analysis comparing:
% DAVID UPDATE -- substantive changes from the original:
%   * Storey + correlation-threshold testing for the observed data is
%     run ONLY over true int-pyr pairs.
%   * Every label permutation runs its OWN Storey + correlation-
%     threshold test over only that draw's pseudo-int x pseudo-pyr pairs.
%   * No global all-vs-all significance mask is created or reused.
%   * No permutation is subsampled to match the real significant-pair count.
%   * The number of significant pseudo pairs is saved for every draw and
%     plotted against the observed number for each animal.

%   H0    = no additional lag shift
%   H+50  = add +50 ms to every significant, oriented lag
%   H-50  = subtract 50 ms from every significant, oriented lag

% UPDATED SIGNIFICANCE/PERMUTATION LOGIC
% --------------------------------------
% REAL DATA:
%   1) Form ONLY the true interneuron x pyramidal-neuron pairs.
%   2) Calculate pair-specific empirical correlation p-values using each
%      pair's circular-shift null correlations.
%   3) Run Storey/mafdr ONLY across those real int-pyr p-values.
%   4) Keep pairs satisfying q <= alpha AND real peak correlation > corrThresh.
%   5) Calculate the real lag imbalance using ALL significant real int-pyr
%      lags. Exact zero lags are excluded from the denominator.

% EACH LABEL-PERMUTATION DRAW:
%   1) Randomly assign nInt neurons as pseudo-interneurons and the remaining
%      nPyr neurons as pseudo-pyramidal neurons.
%   2) Form ONLY the pseudo-int x pseudo-pyr pairs.
%   3) Recalculate pair-specific p-values for that pseudo population.
%   4) Run a NEW Storey/mafdr correction across ONLY that draw's
%      pseudo-int x pseudo-pyr p-values.
%   5) Keep ALL pairs satisfying q <= alpha AND peak correlation > corrThresh.
%   6) Orient each significant lag as pseudo-int -> pseudo-pyr.
%   7) Save the number of significant pairs; DO NOT subsample to match the
%      real significant-pair count.
%   8) Apply the model shift and calculate one lag-imbalance value.

% Thus there is NO universal significance mask and NO fixed-pair-count
% resampling. The real population and every pseudo population receive their
% own Storey + correlation-threshold significance test.

% INPUT:
%   /home/asa7288/pairwise_nochunk_allPairs_ALL_SESSIONS_COMBINED.mat

% OUTPUT:
%   /home/asa7288/pairwiseNoChunkLagImbalanceBayes50ms_StoreyCorrThresh_DAVID_UPDATE.mat

% RUN:
%   computePairwiseNoChunkLagImbalanceBayes50ms_StoreyCorrThresh_DAVID_UPDATE
% or
%   computePairwiseNoChunkLagImbalanceBayes50ms_StoreyCorrThresh_DAVID_UPDATE(0.05, 0.05, 1000)

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
    error(['mafdr was not found. The Bioinformatics Toolbox must be ' ...
        'available on the MATLAB path.']);
end

%% ---------------- settings ----------------

combinedFile = '/home/asa7288/pairwise_nochunk_allPairs_ALL_SESSIONS_COMBINED.mat';
outFile = '/home/asa7288/pairwiseNoChunkLagImbalanceBayes50ms_StoreyCorrThresh_DAVID_UPDATE.mat';
plotDir = '/home/asa7288/pairwiseNoChunkLagImbalance_DAVID_UPDATE_plots';

lagShiftSec = 0.050;

if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

% Reproducible overall run. The H0, H+50, and H-50 helper calls use
% independent draws because the random-number stream continues advancing.
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
results.plotDir = plotDir;
results.alpha = alpha;
results.corrThresh = corrThresh;
results.nNullDraws = nNullDraws;
results.lagShiftSec = lagShiftSec;
results.fdrMethod = [ ...
    'Storey mafdr is run separately on only the relevant cross-population ' ...
    'pair p-values: real int-pyr pairs for the observed statistic and ' ...
    'pseudo-int-pseudo-pyr pairs independently within every permutation draw.'];
results.significanceRule = 'q <= alpha and real peak correlation > corrThresh';
results.lagImbalanceDefinition = [ ...
    '(nPositive - nNegative) / (nPositive + nNegative); exact zero lags excluded'];
results.permutationMethod = [ ...
    'Randomly reassign neuron-type labels while preserving nInt and nPyr; ' ...
    'form all pseudo-int x pseudo-pyr pairs; independently recalculate ' ...
    'pair-specific empirical p-values and Storey q-values for that draw; ' ...
    'retain all significant pairs; orient lags as pseudo-int to pseudo-pyr; ' ...
    'do not match or subsample to the real significant-pair count.'];
results.modelDefinitions = struct( ...
    'H0_shiftSec', 0, ...
    'H50_shiftSec', +lagShiftSec, ...
    'Hneg50_shiftSec', -lagShiftSec);
results.tailProbabilityMethod = [ ...
    'Directional one-sided empirical probabilities with +1 correction. ' ...
    'H0 versus H+50 uses left tails; H0 versus H-50 uses right tails.'];
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

    requiredFields = { ...
        'peakCorrMatAll', 'peakLagMatAll', 'nullCorrMatAllShifts', ...
        'nInt', 'nPyr', 'nAll'};

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
        warning('%s has nInt + nPyr ~= nAll. Skipping session.', animalID);
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
        warning('%s nullCorrMatAllShifts dimensions do not match nAll.', animalID);
        continue;
    end

    fprintf('nInt=%d | nPyr=%d | nAll=%d | shift nulls=%d\n', ...
        nInt, nPyr, nAll, size(nullCorrMat,3));

    %% ====================================================================
    %% real statistic: significance ONLY across true int-pyr pairs
    %% ====================================================================

    % Saved matrix ordering:
    %   1:nInt       = true interneurons
    %   nInt+1:nAll  = true pyramidal neurons
    % Therefore every true int-pyr pair is already stored in the upper
    % triangle and already oriented as real int -> real pyr.
    [actualRows, actualCols] = getIntPyrPairs(nInt, nPyr);

    actualPairStats = computePairSignificance_PairSpecificNull( ...
        peakCorrMat, peakLagMat, nullCorrMat, ...
        actualRows, actualCols, alpha, corrThresh);

    actualSigMask = actualPairStats.sigFDRMask;
    actualSigRows = actualRows(actualSigMask);
    actualSigCols = actualCols(actualSigMask);

    actualSigLags = peakLagMat(sub2ind( ...
        [nAll nAll], actualSigRows, actualSigCols));
    actualSigLags = actualSigLags(:);
    actualSigLags = actualSigLags(isfinite(actualSigLags));

    nActualSigIntPyr = numel(actualSigLags);
    actualLagCounts = countLagSigns(actualSigLags);
    actualLagImbalance = computeLagImbalance(actualSigLags);

    fprintf(['REAL int-pyr tests only: significant=%d/%d | valid tests=%d | ' ...
        'negative=%d | zero=%d | positive=%d | lag imbalance=%.6f\n'], ...
        nActualSigIntPyr, numel(actualRows), actualPairStats.nValidTests, ...
        actualLagCounts.nNegative, actualLagCounts.nZero, ...
        actualLagCounts.nPositive, actualLagImbalance);

    if nActualSigIntPyr == 0
        warning(['%s has no significant real int-pyr pairs. The real lag ' ...
            'imbalance is undefined, but permutation significance-count ' ...
            'distributions will still be generated.'], animalID);
    elseif ~isfinite(actualLagImbalance)
        warning(['%s real lag imbalance is undefined because all significant ' ...
            'real int-pyr lags were exactly zero.'], animalID);
    end

    %% ====================================================================
    %% independently construct H0, H+50, and H-50 distributions
    %% ====================================================================

    fprintf('\nbuilding H0 permutation distribution...\n');
    H0 = computeLabelPermutationLagImbalanceDistribution( ...
        peakCorrMat, peakLagMat, nullCorrMat, ...
        nInt, nPyr, nNullDraws, 0, alpha, corrThresh);

    fprintf('\nbuilding H+50 permutation distribution...\n');
    H50 = computeLabelPermutationLagImbalanceDistribution( ...
        peakCorrMat, peakLagMat, nullCorrMat, ...
        nInt, nPyr, nNullDraws, +lagShiftSec, alpha, corrThresh);

    fprintf('\nbuilding H-50 permutation distribution...\n');
    Hneg50 = computeLabelPermutationLagImbalanceDistribution( ...
        peakCorrMat, peakLagMat, nullCorrMat, ...
        nInt, nPyr, nNullDraws, -lagShiftSec, alpha, corrThresh);

    %% ---------- retain valid lag-imbalance null values ----------

    validH0 = H0.lagImbalance(isfinite(H0.lagImbalance));
    validH50 = H50.lagImbalance(isfinite(H50.lagImbalance));
    validHneg50 = Hneg50.lagImbalance(isfinite(Hneg50.lagImbalance));

    %% ====================================================================
    %% directional one-sided empirical probabilities
    %% ====================================================================

    if isfinite(actualLagImbalance)
        % H+50 predicts positive lag imbalance; use left tails.
        pValH0_left = empiricalLeftTailProbability(actualLagImbalance, validH0);
        pValH50_left = empiricalLeftTailProbability(actualLagImbalance, validH50);

        % H-50 predicts negative lag imbalance; use right tails.
        pValH0_right = empiricalRightTailProbability(actualLagImbalance, validH0);
        pValHneg50_right = empiricalRightTailProbability(actualLagImbalance, validHneg50);
    else
        pValH0_left = NaN;
        pValH50_left = NaN;
        pValH0_right = NaN;
        pValHneg50_right = NaN;
    end

    evidenceRatio_H0_over_H50 = safeRatio(pValH0_left, pValH50_left);
    evidenceRatio_H0_over_Hneg50 = safeRatio(pValH0_right, pValHneg50_right);

    %% ---------- descriptive intervals ----------

    nullCIH0 = safePercentileInterval(validH0, [2.5 97.5]);
    nullCIH50 = safePercentileInterval(validH50, [2.5 97.5]);
    nullCIHneg50 = safePercentileInterval(validHneg50, [2.5 97.5]);

    sigCountCIH0 = safePercentileInterval(H0.nSignificantPairs, [2.5 97.5]);
    sigCountCIH50 = safePercentileInterval(H50.nSignificantPairs, [2.5 97.5]);
    sigCountCIHneg50 = safePercentileInterval(Hneg50.nSignificantPairs, [2.5 97.5]);

    %% ====================================================================
    %% plot number-significant-pairs and lag-imbalance distributions
    %% ====================================================================

    plotBaseName = sprintf('%s_UPDATED_significance_and_lag_imbalance', animalID);
    pngFile = fullfile(plotDir, [plotBaseName '.png']);
    figFile = fullfile(plotDir, [plotBaseName '.fig']);

    createSessionDiagnosticFigure( ...
        animalID, nActualSigIntPyr, actualLagImbalance, ...
        H0, H50, Hneg50, pngFile, figFile);

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
    R.lagShiftSec = lagShiftSec;

    % Real int-pyr significance/statistic
    R.actualRows = actualRows;
    R.actualCols = actualCols;
    R.actualPairStats = actualPairStats;
    R.actualSigMask = actualSigMask;
    R.actualSigRows = actualSigRows;
    R.actualSigCols = actualSigCols;
    R.actualSigLags = actualSigLags;
    R.nActualSigIntPyr = nActualSigIntPyr;
    R.actualLagCounts = actualLagCounts;
    R.actualLagImbalance = actualLagImbalance;

    % Full per-draw permutation results, including significance counts
    R.H0 = H0;
    R.H50 = H50;
    R.Hneg50 = Hneg50;

    % Convenience copies
    R.nullLagImbalanceH0 = H0.lagImbalance;
    R.nullLagImbalanceH50 = H50.lagImbalance;
    R.nullLagImbalanceHneg50 = Hneg50.lagImbalance;

    R.validNullLagImbalanceH0 = validH0;
    R.validNullLagImbalanceH50 = validH50;
    R.validNullLagImbalanceHneg50 = validHneg50;

    R.nSignificantPairsH0 = H0.nSignificantPairs;
    R.nSignificantPairsH50 = H50.nSignificantPairs;
    R.nSignificantPairsHneg50 = Hneg50.nSignificantPairs;

    R.nullCIH0 = nullCIH0;
    R.nullCIH50 = nullCIH50;
    R.nullCIHneg50 = nullCIHneg50;

    R.sigCountCIH0 = sigCountCIH0;
    R.sigCountCIH50 = sigCountCIH50;
    R.sigCountCIHneg50 = sigCountCIHneg50;

    R.pValH0_left = pValH0_left;
    R.pValH50_left = pValH50_left;
    R.pValH0_right = pValH0_right;
    R.pValHneg50_right = pValHneg50_right;

    R.evidenceRatio_H0_over_H50 = evidenceRatio_H0_over_H50;
    R.evidenceRatio_H0_over_Hneg50 = evidenceRatio_H0_over_Hneg50;

    R.plotPNG = pngFile;
    R.plotFIG = figFile;

    R.modelConstruction = [ ...
        'Real significance is calculated only across true int-pyr pairs. ' ...
        'Each null draw independently reassigns pseudo labels, forms only ' ...
        'pseudo-int x pseudo-pyr pairs, recalculates pair-specific p-values, ' ...
        'runs Storey mafdr only across those p-values, applies the correlation ' ...
        'threshold, retains every significant pair, orients its lag as ' ...
        'pseudo-int to pseudo-pyr, applies the model shift, and calculates ' ...
        'lag imbalance. No universal significance mask or pair-count matching ' ...
        'is used.'];

    results.sessions{sessInd} = R;

    %% ---------- print session summary ----------

    fprintf('\n%s UPDATED summary:\n', animalID);
    fprintf('  real significant int-pyr pairs: %d/%d\n', ...
        nActualSigIntPyr, numel(actualRows));
    fprintf('  real lag imbalance: %.6f\n', actualLagImbalance);
    fprintf(['  permutation significant-pair counts (median [2.5,97.5]):\n' ...
        '    H0:    %.1f [%.1f, %.1f]\n' ...
        '    H+50:  %.1f [%.1f, %.1f]\n' ...
        '    H-50:  %.1f [%.1f, %.1f]\n'], ...
        median(H0.nSignificantPairs), sigCountCIH0(1), sigCountCIH0(2), ...
        median(H50.nSignificantPairs), sigCountCIH50(1), sigCountCIH50(2), ...
        median(Hneg50.nSignificantPairs), sigCountCIHneg50(1), sigCountCIHneg50(2));
    fprintf('  valid lag-imbalance draws: H0=%d | H+50=%d | H-50=%d\n', ...
        numel(validH0), numel(validH50), numel(validHneg50));
    fprintf('  evidence ratios: H0/H+50=%.6f | H0/H-50=%.6f\n', ...
        evidenceRatio_H0_over_H50, evidenceRatio_H0_over_Hneg50);
    fprintf('  saved diagnostic plot: %s\n', pngFile);
end

%% ---------------- save all sessions ----------------

save(outFile, 'results', '-v7.3');

fprintf('\n========================================\n');
fprintf('saved:\n%s\n', outFile);
fprintf('plots saved in:\n%s\n', plotDir);
fprintf('========================================\n');

end

%% ========================================================================
%% helper: true int-pyr pair indices
%% ========================================================================

function [rows, cols] = getIntPyrPairs(nInt, nPyr)
% Saved neuron ordering:
%   1:nInt             = true interneurons
%   nInt+1:nInt+nPyr   = true pyramidal neurons
% Therefore all true int-pyr pairs are in the upper triangle and already
% have the desired int -> pyr orientation.

[rr, cc] = ndgrid(1:nInt, nInt + (1:nPyr));
rows = rr(:);
cols = cc(:);

end

%% ========================================================================
%% helper: calculate Storey + correlation-threshold significance for an
%% arbitrary supplied list of upper-triangle neuron pairs
%% ========================================================================

function out = computePairSignificance_PairSpecificNull( ...
    realMat, lagMat, nullMat, rows, cols, alpha, corrThresh)

nPairs = numel(rows);
realVals = nan(nPairs,1);
lagVals = nan(nPairs,1);
pVals = nan(nPairs,1);
nNullPerPair = zeros(nPairs,1);

for pairInd = 1:nPairs
    r = rows(pairInd);
    c = cols(pairInd);

    if r >= c
        error('Pair indices supplied to significance helper must satisfy row < column.');
    end

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
    % How often is the circular-shift null correlation at least as large as
    % the observed real peak correlation for this same neuron pair?
    pVals(pairInd) = ...
        (sum(thisNull >= thisReal) + 1) / (numel(thisNull) + 1);
end

validP = isfinite(pVals);
qVals = nan(size(pVals));
pFDRVals = nan(size(pVals));

if any(validP)
    [pFDRtmp, qTmp] = mafdr(pVals(validP));
    pFDRVals(validP) = pFDRtmp;
    qVals(validP) = qTmp;
end

sigFDRMask = false(size(pVals));
sigFDRMask(validP) = ...
    (qVals(validP) <= alpha) & (realVals(validP) > corrThresh);

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
out.nSigFDR = nnz(sigFDRMask);
out.sigFDRMask = sigFDRMask;

end

%% ========================================================================
%% helper: independently significance-test every pseudo population and build
%% one model-specific lag-imbalance distribution
%% ========================================================================

function out = computeLabelPermutationLagImbalanceDistribution( ...
    peakCorrMat, peakLagMat, nullCorrMat, ...
    nInt, nPyr, nNullDraws, modelShiftSec, alpha, corrThresh)

nAll = nInt + nPyr;
maxPossiblePairs = nInt * nPyr;

out = struct();
out.modelShiftSec = modelShiftSec;
out.lagImbalance = nan(nNullDraws,1);
out.nSignificantPairs = zeros(nNullDraws,1);
out.nValidTests = zeros(nNullDraws,1);
out.nPairsTested = repmat(maxPossiblePairs, nNullDraws, 1);
out.selectedPositiveCounts = zeros(nNullDraws,1);
out.selectedNegativeCounts = zeros(nNullDraws,1);
out.selectedZeroCounts = zeros(nNullDraws,1);

for drawInd = 1:nNullDraws

    %% ---------- randomly reassign neuron-type labels ----------

    permutedNeuronOrder = randperm(nAll);
    pseudoIntInds = permutedNeuronOrder(1:nInt);
    pseudoPyrInds = permutedNeuronOrder(nInt+1:end);

    %% ---------- form every pseudo-int x pseudo-pyr pair ----------

    [pseudoIntGrid, pseudoPyrGrid] = ndgrid(pseudoIntInds, pseudoPyrInds);
    pseudoIntForPair = pseudoIntGrid(:);
    pseudoPyrForPair = pseudoPyrGrid(:);

    % Convert each conceptual pseudo-int -> pseudo-pyr pair into the saved
    % upper-triangle matrix location. orientationSign converts the stored
    % row->column lag back into pseudo-int->pseudo-pyr orientation.
    upperRows = min(pseudoIntForPair, pseudoPyrForPair);
    upperCols = max(pseudoIntForPair, pseudoPyrForPair);

    orientationSigns = ones(maxPossiblePairs,1);
    orientationSigns(pseudoIntForPair > pseudoPyrForPair) = -1;

    %% ---------- NEW: significance test for this pseudo population only ----------

    pairStats = computePairSignificance_PairSpecificNull( ...
        peakCorrMat, peakLagMat, nullCorrMat, ...
        upperRows, upperCols, alpha, corrThresh);

    sigMask = pairStats.sigFDRMask;

    sigUpperRows = upperRows(sigMask);
    sigUpperCols = upperCols(sigMask);
    sigOrientationSigns = orientationSigns(sigMask);

    storedSigLags = peakLagMat(sub2ind( ...
        [nAll nAll], sigUpperRows, sigUpperCols));

    orientedSigLags = sigOrientationSigns .* storedSigLags;
    finiteMask = isfinite(orientedSigLags);

    sigUpperRows = sigUpperRows(finiteMask);
    sigUpperCols = sigUpperCols(finiteMask);
    sigOrientationSigns = sigOrientationSigns(finiteMask);
    orientedSigLags = orientedSigLags(finiteMask);

    % Use every significant pseudo-int-pseudo-pyr pair. There is no
    % subsampling and no attempt to match the real significant-pair count.
    modelLags = orientedSigLags + modelShiftSec;
    thisLagImbalance = computeLagImbalance(modelLags);
    thisCounts = countLagSigns(modelLags);

    %% ---------- save this draw ----------

    out.lagImbalance(drawInd) = thisLagImbalance;
    out.nSignificantPairs(drawInd) = numel(orientedSigLags);
    out.nValidTests(drawInd) = pairStats.nValidTests;
    out.selectedPositiveCounts(drawInd) = thisCounts.nPositive;
    out.selectedNegativeCounts(drawInd) = thisCounts.nNegative;
    out.selectedZeroCounts(drawInd) = thisCounts.nZero;


    if mod(drawInd,10) == 0 || drawInd == nNullDraws
        fprintf(['  completed draw %d/%d | model shift=%+.3f s | ' ...
            'significant pairs=%d/%d | lag imbalance=%.6f\n'], ...
            drawInd, nNullDraws, modelShiftSec, ...
            out.nSignificantPairs(drawInd), maxPossiblePairs, ...
            thisLagImbalance);
    end
end

end

%% ========================================================================
%% helper: diagnostic figure for each animal
%% ========================================================================

function createSessionDiagnosticFigure( ...
    animalID, nActualSigPairs, actualLagImbalance, ...
    H0, H50, Hneg50, pngFile, figFile)

fig = figure('Color', 'w', 'Visible', 'off', ...
    'Name', sprintf('%s updated pairwise permutation diagnostics', animalID), ...
    'Position', [100 100 1350 520]);

T = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(T, sprintf('%s: per-permutation Storey testing on pseudo-int x pseudo-pyr pairs', animalID), ...
    'Interpreter', 'none');

%% Significant-pair count distribution
nexttile;
hold on;

allCounts = [H0.nSignificantPairs(:); H50.nSignificantPairs(:); Hneg50.nSignificantPairs(:)];
countEdges = integerHistogramEdges(allCounts);

histogram(H0.nSignificantPairs, countEdges, ...
    'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 1.6);
histogram(H50.nSignificantPairs, countEdges, ...
    'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 1.6);
histogram(Hneg50.nSignificantPairs, countEdges, ...
    'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 1.6);
xline(nActualSigPairs, '--', 'LineWidth', 2);

xlabel('Number of significant pseudo-int-pseudo-pyr pairs');
ylabel('Probability');
title(sprintf('Significant-pair counts | real = %d', nActualSigPairs));
legend({'H0 permutations', 'H+50 permutations', 'H-50 permutations', ...
    'Real significant int-pyr count'}, 'Location', 'best');
grid on;
box off;

%% Lag-imbalance null distributions
nexttile;
hold on;

lagEdges = linspace(-1, 1, 31);
validH0 = H0.lagImbalance(isfinite(H0.lagImbalance));
validH50 = H50.lagImbalance(isfinite(H50.lagImbalance));
validHneg50 = Hneg50.lagImbalance(isfinite(Hneg50.lagImbalance));

histogram(validH0, lagEdges, ...
    'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 1.6);
histogram(validH50, lagEdges, ...
    'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 1.6);
histogram(validHneg50, lagEdges, ...
    'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 1.6);

if isfinite(actualLagImbalance)
    xline(actualLagImbalance, '--', 'LineWidth', 2);
end

xlabel('Lag imbalance');
ylabel('Probability');
title(sprintf('Null lag-imbalance distributions | real = %.3f', actualLagImbalance));
legend({'H0', 'H+50', 'H-50', 'Real lag imbalance'}, 'Location', 'best');
xlim([-1 1]);
grid on;
box off;

try
    exportgraphics(fig, pngFile, 'Resolution', 300);
catch
    saveas(fig, pngFile);
end

savefig(fig, figFile);
close(fig);

end

%% ========================================================================
%% helper: integer-centered histogram edges
%% ========================================================================

function edges = integerHistogramEdges(values)

values = values(isfinite(values));

if isempty(values)
    edges = [-0.5 0.5];
    return;
end

minVal = floor(min(values));
maxVal = ceil(max(values));

if minVal == maxVal
    edges = [minVal - 0.5, maxVal + 0.5];
else
    edges = (minVal - 0.5):1:(maxVal + 0.5);
end

end

%% ========================================================================
%% helper: lag imbalance
%% ========================================================================

function lagImbalance = computeLagImbalance(lags)
% Exact zero lags are excluded from the denominator.
% Range:
%   -1 = all nonzero lags are negative
%    0 = equal positive and negative counts
%   +1 = all nonzero lags are positive

lags = lags(isfinite(lags));
nPositive = nnz(lags > 0);
nNegative = nnz(lags < 0);
denominator = nPositive + nNegative;

if denominator == 0
    lagImbalance = NaN;
else
    lagImbalance = (nPositive - nNegative) / denominator;
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
%% helper: empirical left-tail probability
%% ========================================================================

function pVal = empiricalLeftTailProbability(actualValue, nullValues)

nullValues = nullValues(isfinite(nullValues));

if ~isfinite(actualValue) || isempty(nullValues)
    pVal = NaN;
    return;
end

pVal = (sum(nullValues <= actualValue) + 1) / (numel(nullValues) + 1);

end

%% ========================================================================
%% helper: empirical right-tail probability
%% ========================================================================

function pVal = empiricalRightTailProbability(actualValue, nullValues)

nullValues = nullValues(isfinite(nullValues));

if ~isfinite(actualValue) || isempty(nullValues)
    pVal = NaN;
    return;
end

pVal = (sum(nullValues >= actualValue) + 1) / (numel(nullValues) + 1);

end

%% ========================================================================
%% helper: safe evidence-ratio calculation
%% ========================================================================

function ratio = safeRatio(numerator, denominator)

if ~isfinite(numerator) || ~isfinite(denominator) || denominator == 0
    ratio = NaN;
else
    ratio = numerator / denominator;
end

end

%% ========================================================================
%% helper: percentile interval that tolerates empty input
%% ========================================================================

function interval = safePercentileInterval(values, requestedPercentiles)

values = values(isfinite(values));

if isempty(values)
    interval = [NaN NaN];
else
    interval = prctile(values, requestedPercentiles);
end

end
