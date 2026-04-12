function runIntPyrLagImbalanceMAFDRPermutationTest(saveFile, analysisType, alpha, nNullDraws, corrThresh)
% runs int-vs-pyr-only storey/mafdr analysis and lag-imbalance null-interval testing for either no-chunk or chunked pairwise saved results

% updated logic:
%   - compute empirical p-values for all true int x pyr pairs
%   - run matlab mafdr on all valid p-values
%   - call pairs significant if: q <= alpha and peak correlation > corrThresh
%   - compute lag imbalance from q-significant lags: (# positive lags - # negative lags) / (# positive lags + # negative lags) , zero lags are ignored

% null test:
%   - sample the same number of unique pairs from all-vs-all upper-triangle
%   - sampled pairs may include int-int, pyr-pyr, and int-pyr
%   - no self-pairs and no duplicate pairs within a draw
%   - for null draws only, randomly assign each sampled pair a sign of +1 or -1 before using its lag, to account for the fact that reversing pair ordercwould flip the lag sign
%   - compute q-significant lags and lag imbalance
%   - repeat nNullDraws times
%   - determine significance directly from whether actual lag imbalance fallscoutside the central 95% null interval

% inputs
%   saveFile : combined pairwise .mat file
%   analysisType : "nochunk" or "chunked"
%   alpha : q-value threshold (default 0.05)
%   nNullDraws : number of null lag-imbalance draws (default 100)
%   corrThresh : minimum peak correlation threshold (default 0.05)

if nargin < 3 || isempty(alpha)
    alpha = 0.05;
end
if nargin < 4 || isempty(nNullDraws)
    nNullDraws = 100;
end
if nargin < 5 || isempty(corrThresh)
    corrThresh = 0.05;
end

analysisType = lower(string(analysisType));

if exist('mafdr', 'file') ~= 2
    error('mafdr is not on the matlab path / bioinformatics toolbox is unavailable.');
end

S = load(char(saveFile));

switch analysisType
    case "nochunk"
        if ~isfield(S, 'allSessions')
            error('nochunk file must contain allSessions.');
        end
        nSess = numel(S.allSessions.sessions);

    case "chunked"
        if ~isfield(S, 'all_peakCorrMat_all')
            error('chunked file structure not recognized.');
        end
        nSess = numel(S.all_peakCorrMat_all);

    otherwise
        error('analysisType must be "nochunk" or "chunked".');
end

rng(0);

results = struct();
results.sourceFile = saveFile;
results.analysisType = char(analysisType);
results.alpha = alpha;
results.nNullDraws = nNullDraws;
results.corrThresh = corrThresh;
results.fdrMethod = 'mafdr_storey_all_valid_pvalues_plus_corr_threshold';
results.directionMetric = 'lagImbalance=(nPositive-nNegative)/(nPositive+nNegative), zero lags ignored';
results.significanceRule = 'actual lag imbalance outside central 95% null interval';
results.nullLagSignHandling = 'null sampled upper-triangle pairs assigned random +/- sign to emulate both pair orientations';
results.sessions = cell(1, nSess);

for s = 1:nSess

    [peakCorr, peakLag, nullCorr, nInt, nPyr, nAll, animalID] = getSessionData(S, s, analysisType);

    fprintf('\n=============================\n');
    fprintf('processing %s (%s)\n', animalID, analysisType);
    fprintf('nInt = %d | nPyr = %d | nAll = %d\n', nInt, nPyr, nAll);
    fprintf('=============================\n');

    % ---------------- actual int x pyr pairs ----------------
    [actualRows, actualCols] = getIntPyrPairs(nInt, nPyr);

    % actual keeps its consistent int->pyr sign convention
    actual = computePairSetStatsMAFDR(peakCorr, peakLag, nullCorr, actualRows, actualCols, alpha, corrThresh);

    fprintf('%s actual: nominal pairs = %d | valid tests = %d | uncorr p<=%.3f = %d | q<=%.3f & corr>%.3f sig = %d | lag imbalance = %.6f\n', ...
        animalID, numel(actualRows), actual.nValidTests, alpha, actual.nSigUncorr, alpha, corrThresh, actual.nSigFDR, actual.lagImbalance);

    % ---------------- null draws from all-vs-all upper triangle ----------------
    [poolRows, poolCols] = getAllUpperPairs(nAll);

    nActualNominalPairs = numel(actualRows);

    if nActualNominalPairs > numel(poolRows)
        error('%s: actual nominal pair count exceeds all-vs-all upper-triangle pool size.', animalID);
    end

    nullLagImbalance = nan(nNullDraws,1);
    nullNSigUncorr = nan(nNullDraws,1);
    nullNSigFDR = nan(nNullDraws,1);
    nullNValidTests = nan(nNullDraws,1);
    nullLagCell = cell(nNullDraws,1);
    nullPairTypeCounts = nan(nNullDraws, 3); % [int-int, int-pyr, pyr-pyr]
    nullLagSignVecCell = cell(nNullDraws,1);

    for r = 1:nNullDraws
        drawIdx = randperm(numel(poolRows), nActualNominalPairs);
        drawRows = poolRows(drawIdx);
        drawCols = poolCols(drawIdx);

        % key fix:
        % for null sampled pairs, randomly assign an orientation sign
        % because using only upper-triangle ordering otherwise biases lag sign
        lagSignVec = ones(nActualNominalPairs,1);
        lagSignVec(rand(nActualNominalPairs,1) > 0.5) = -1;

        nullDraw = computePairSetStatsMAFDR( ...
            peakCorr, peakLag, nullCorr, drawRows, drawCols, alpha, corrThresh, lagSignVec);

        nullLagImbalance(r) = nullDraw.lagImbalance;
        nullNSigUncorr(r) = nullDraw.nSigUncorr;
        nullNSigFDR(r) = nullDraw.nSigFDR;
        nullNValidTests(r) = nullDraw.nValidTests;
        nullLagCell{r} = nullDraw.sigLagVec;
        nullLagSignVecCell{r} = lagSignVec;

        nullPairTypeCounts(r,:) = classifySampledPairs(drawRows, drawCols, nInt);
    end

    validNullLagImbalance = nullLagImbalance(~isnan(nullLagImbalance) & isfinite(nullLagImbalance));

    if ~isempty(validNullLagImbalance)
        lagImbalanceNullCI = prctile(validNullLagImbalance, [2.5 97.5]);
    else
        lagImbalanceNullCI = [NaN NaN];
    end

    lagImbalanceSigByCI = false;
    if ~isnan(actual.lagImbalance) && all(~isnan(lagImbalanceNullCI))
        lagImbalanceSigByCI = ...
            (actual.lagImbalance < lagImbalanceNullCI(1)) || ...
            (actual.lagImbalance > lagImbalanceNullCI(2));
    end

    fprintf('%s lag imbalance test: actual = %.6f | null n = %d | 95%% CI = [%.6f, %.6f] | sig by CI = %d\n', ...
        animalID, actual.lagImbalance, numel(validNullLagImbalance), ...
        lagImbalanceNullCI(1), lagImbalanceNullCI(2), lagImbalanceSigByCI);

    % ---------------- plots ----------------

    figure('Name', sprintf('%s actual int-pyr significant lags (mafdr)', animalID), 'Color', 'w');
    if ~isempty(actual.sigLagVec)
        lagEdgesActual = makeLagEdges(actual.sigLagVec);
        histogram(actual.sigLagVec, 'BinEdges', lagEdgesActual, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    else
        histogram(nan, 'BinEdges', -0.201:0.001:0.201, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    end
    xlabel('Peak Lag (s)');
    ylabel('Count');
    title(sprintf('%s actual int-pyr q<=%.3f & corr>%.3f lags | uncorr=%d | sig=%d | lag imbalance=%.4f', ...
        animalID, alpha, corrThresh, actual.nSigUncorr, actual.nSigFDR, actual.lagImbalance));
    grid on;

    pooledNullLags = vertcat(nullLagCell{:});
    pooledNullLags = pooledNullLags(~isnan(pooledNullLags) & isfinite(pooledNullLags));

    figure('Name', sprintf('%s null significant lags (mafdr)', animalID), 'Color', 'w');
    if ~isempty(pooledNullLags)
        lagEdgesNull = makeLagEdges(pooledNullLags);
        histogram(pooledNullLags, 'BinEdges', lagEdgesNull, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    else
        histogram(nan, 'BinEdges', -0.201:0.001:0.201, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    end
    xlabel('Peak Lag (s)');
    ylabel('Count');
    title(sprintf('%s pooled null significant lags across %d draws', animalID, nNullDraws));
    grid on;

    figure('Name', sprintf('%s lag imbalance null distribution (mafdr)', animalID), 'Color', 'w');
    if ~isempty(validNullLagImbalance)
        histogram(validNullLagImbalance, 30, 'FaceAlpha', 0.8, 'EdgeColor', 'none'); hold on;
        xline(actual.lagImbalance, 'r-', 'LineWidth', 2);
        xline(lagImbalanceNullCI(1), 'k--', 'LineWidth', 1.5);
        xline(lagImbalanceNullCI(2), 'k--', 'LineWidth', 1.5);
        legend({'Null lag imbalance', 'Actual lag imbalance', '95% null interval'}, 'Location', 'best');
    else
        histogram(nan, 30, 'FaceAlpha', 0.8, 'EdgeColor', 'none'); hold on;
        xline(actual.lagImbalance, 'r-', 'LineWidth', 2);
        legend({'Null lag imbalance', 'Actual lag imbalance'}, 'Location', 'best');
    end
    xlabel('Lag Imbalance = (nPositive - nNegative) / (nPositive + nNegative)');
    ylabel('Count');
    title(sprintf('%s lag imbalance null distribution | actual=%.4f | sig by 95%% CI = %d', ...
        animalID, actual.lagImbalance, lagImbalanceSigByCI));
    grid on;

    % ---------------- save session results ----------------
    R = struct();
    R.animalID = animalID;
    R.nInt = nInt;
    R.nPyr = nPyr;
    R.nAll = nAll;

    R.actualRows = actualRows;
    R.actualCols = actualCols;
    R.actual = actual;
    R.actual.nPairsNominal = nActualNominalPairs;

    R.nullLagImbalance = nullLagImbalance;
    R.validNullLagImbalance = validNullLagImbalance;
    R.lagImbalanceNullCI = lagImbalanceNullCI;
    R.lagImbalanceSigByCI = lagImbalanceSigByCI;

    R.nullNSigUncorr = nullNSigUncorr;
    R.nullNSigFDR = nullNSigFDR;
    R.nullNValidTests = nullNValidTests;
    R.nullLagCell = nullLagCell;
    R.nullPairTypeCounts = nullPairTypeCounts;
    R.nullLagSignVecCell = nullLagSignVecCell;

    results.sessions{s} = R;
end

[folderPath, baseName, ~] = fileparts(char(saveFile));
outFile = fullfile(folderPath, sprintf('%s_intPyrLagImbalanceMAFDR_%s.mat', baseName, char(analysisType)));
save(char(outFile), 'results', '-v7.3');

fprintf('\nsaved:\n%s\n', outFile);

end

% ============================================================
function [peakCorr, peakLag, nullCorr, nInt, nPyr, nAll, animalID] = getSessionData(S, s, analysisType)

switch analysisType
    case "nochunk"
        sess = S.allSessions.sessions(s);

        peakCorr = sess.peakCorrMatAll;
        peakLag = sess.peakLagMatAll;
        nullCorr = sess.nullCorrMatAllShifts;
        nInt = sess.nInt;
        nPyr = sess.nPyr;
        nAll = sess.nAll;

        animalID = regexp(sess.baseDir, 'D\d+', 'match', 'once');
        if isempty(animalID)
            animalID = sprintf('Session%d', s);
        end

    case "chunked"
        peakCorr = S.all_peakCorrMat_all{s};
        peakLag = S.all_peakLagSecMat_all{s};
        nullCorr = S.all_nullCorrMat_allShifts{s};
        nInt = S.all_nInt_ref(s);
        nPyr = S.all_nPyr_ref(s);
        nAll = S.all_nAll(s);

        animalID = regexp(S.baseDirs{s}, 'D\d+', 'match', 'once');
        if isempty(animalID)
            animalID = sprintf('Session%d', s);
        end
end
end

% ============================================================
function [rows, cols] = getIntPyrPairs(nInt, nPyr)
[rr, cc] = ndgrid(1:nInt, nInt + (1:nPyr));
rows = rr(:);
cols = cc(:);
end

% ============================================================
function [rows, cols] = getAllUpperPairs(nAll)
mask = triu(true(nAll), 1);
[rows, cols] = find(mask);
end

% ============================================================
function out = computePairSetStatsMAFDR(realMat, lagMat, nullMat, rows, cols, alpha, corrThresh, lagSignVec)

if nargin < 8 || isempty(lagSignVec)
    lagSignVec = ones(numel(rows),1);
end

nPairs = numel(rows);

if numel(lagSignVec) ~= nPairs
    error('lagSignVec must have one entry per selected pair.');
end

realVals = nan(nPairs,1);
lagVals = nan(nPairs,1);
pVals = nan(nPairs,1);

for i = 1:nPairs
    r = rows(i);
    c = cols(i);

    thisReal = realMat(r,c);
    thisLag = lagMat(r,c);
    thisNull = squeeze(nullMat(r,c,:));
    thisNull = thisNull(~isnan(thisNull) & isfinite(thisNull));

    realVals(i) = thisReal;
    lagVals(i) = lagSignVec(i) * thisLag;

    if isnan(thisReal) || ~isfinite(thisReal) || isnan(thisLag) || ~isfinite(thisLag) || numel(thisNull) < 10
        pVals(i) = NaN;
        continue;
    end

    pVals(i) = (sum(thisNull >= thisReal) + 1) / (numel(thisNull) + 1);
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
out.lagImbalance = computeLagImbalance(sigLagVec);
out.nPositiveLag = sum(sigLagVec > 0);
out.nNegativeLag = sum(sigLagVec < 0);
out.nZeroLag = sum(sigLagVec == 0);
out.sigUncorrMask = sigUncorr;
out.sigFDRMask = sigFDR;
out.lagSignVec = lagSignVec;
end

% ============================================================
function val = computeLagImbalance(x)
x = x(~isnan(x) & isfinite(x));

if isempty(x)
    val = NaN;
    return;
end

nPos = sum(x > 0);
nNeg = sum(x < 0);
den = nPos + nNeg;

if den == 0
    val = NaN;
else
    val = (nPos - nNeg) / den;
end
end

% ============================================================
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

% ============================================================
function counts = classifySampledPairs(rows, cols, nInt)
isIntRow = rows <= nInt;
isIntCol = cols <= nInt;

nIntInt = sum(isIntRow & isIntCol);
nPyrPyr = sum(~isIntRow & ~isIntCol);
nIntPyr = sum(xor(isIntRow, isIntCol));

counts = [nIntInt, nIntPyr, nPyrPyr];
end
