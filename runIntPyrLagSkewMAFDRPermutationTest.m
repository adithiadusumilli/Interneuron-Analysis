function runIntPyrLagSkewMAFDRPermutationTest(saveFile, analysisType, alpha, nNullDraws, corrThresh)
% runs int-vs-pyr-only storey/mafdr analysis and skew permutation testing for either no-chunk or chunked pairwise saved results

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
results.directionMetric = 'skew=(mean-median)/std of significant lag distribution';
results.significanceRule = 'actual skew outside central 95% null interval';
results.nullLagSignHandling = 'null sampled upper-triangle pairs assigned random +/- sign to emulate both pair orientations';
results.sessions = cell(1, nSess);

for s = 1:nSess

    [peakCorr, peakLag, nullCorr, nInt, nPyr, nAll, animalID] = getSessionData(S, s, analysisType);

    fprintf('\n=============================\n');
    fprintf('processing %s (%s)\n', animalID, analysisType);
    fprintf('nInt = %d | nPyr = %d | nAll = %d\n', nInt, nPyr, nAll);
    fprintf('=============================\n');

    [actualRows, actualCols] = getIntPyrPairs(nInt, nPyr);

    actual = computePairSetStatsMAFDRSkew(peakCorr, peakLag, nullCorr, actualRows, actualCols, alpha, corrThresh);

    fprintf('%s actual: nominal pairs = %d | valid tests = %d | uncorr p<=%.3f = %d | q<=%.3f & corr>%.3f sig = %d | skew = %.6f\n', ...
        animalID, numel(actualRows), actual.nValidTests, alpha, actual.nSigUncorr, alpha, corrThresh, actual.nSigFDR, actual.skew);

    [poolRows, poolCols] = getAllUpperPairs(nAll);
    nActualNominalPairs = numel(actualRows);

    if nActualNominalPairs > numel(poolRows)
        error('%s: actual nominal pair count exceeds all-vs-all upper-triangle pool size.', animalID);
    end

    nullSkews = nan(nNullDraws,1);
    nullNSigUncorr = nan(nNullDraws,1);
    nullNSigFDR = nan(nNullDraws,1);
    nullNValidTests = nan(nNullDraws,1);
    nullLagCell = cell(nNullDraws,1);
    nullPairTypeCounts = nan(nNullDraws, 3);
    nullLagSignVecCell = cell(nNullDraws,1);

    for r = 1:nNullDraws
        drawIdx = randperm(numel(poolRows), nActualNominalPairs);
        drawRows = poolRows(drawIdx);
        drawCols = poolCols(drawIdx);

        lagSignVec = ones(nActualNominalPairs,1);
        lagSignVec(rand(nActualNominalPairs,1) > 0.5) = -1;

        nullDraw = computePairSetStatsMAFDRSkew( ...
            peakCorr, peakLag, nullCorr, drawRows, drawCols, alpha, corrThresh, lagSignVec);

        nullSkews(r) = nullDraw.skew;
        nullNSigUncorr(r) = nullDraw.nSigUncorr;
        nullNSigFDR(r) = nullDraw.nSigFDR;
        nullNValidTests(r) = nullDraw.nValidTests;
        nullLagCell{r} = nullDraw.sigLagVec;
        nullLagSignVecCell{r} = lagSignVec;

        nullPairTypeCounts(r,:) = classifySampledPairs(drawRows, drawCols, nInt);
    end

    validNullSkews = nullSkews(~isnan(nullSkews) & isfinite(nullSkews));

    if ~isempty(validNullSkews)
        skewNullCI = prctile(validNullSkews, [2.5 97.5]);
    else
        skewNullCI = [NaN NaN];
    end

    skewSigByCI = false;
    if ~isnan(actual.skew) && all(~isnan(skewNullCI))
        skewSigByCI = (actual.skew < skewNullCI(1)) || (actual.skew > skewNullCI(2));
    end

    fprintf('%s skew test: actual = %.6f | null n = %d | 95%% CI = [%.6f, %.6f] | sig by CI = %d\n', ...
        animalID, actual.skew, numel(validNullSkews), ...
        skewNullCI(1), skewNullCI(2), skewSigByCI);

    R = struct();
    R.animalID = animalID;
    R.nInt = nInt;
    R.nPyr = nPyr;
    R.nAll = nAll;

    R.actualRows = actualRows;
    R.actualCols = actualCols;
    R.actual = actual;
    R.actual.nPairsNominal = nActualNominalPairs;

    R.nullSkews = nullSkews;
    R.validNullSkews = validNullSkews;
    R.skewNullCI = skewNullCI;
    R.skewSigByCI = skewSigByCI;

    R.nullNSigUncorr = nullNSigUncorr;
    R.nullNSigFDR = nullNSigFDR;
    R.nullNValidTests = nullNValidTests;
    R.nullLagCell = nullLagCell;
    R.nullPairTypeCounts = nullPairTypeCounts;
    R.nullLagSignVecCell = nullLagSignVecCell;

    results.sessions{s} = R;
end

[folderPath, baseName, ~] = fileparts(char(saveFile));
outFile = fullfile(folderPath, sprintf('%s_intPyrSkewMAFDR_%s.mat', baseName, char(analysisType)));
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
function out = computePairSetStatsMAFDRSkew(realMat, lagMat, nullMat, rows, cols, alpha, corrThresh, lagSignVec)

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
out.skew = computeSkew(sigLagVec);
out.sigUncorrMask = sigUncorr;
out.sigFDRMask = sigFDR;
out.lagSignVec = lagSignVec;
end

% ============================================================
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

% ============================================================
function counts = classifySampledPairs(rows, cols, nInt)
isIntRow = rows <= nInt;
isIntCol = cols <= nInt;

nIntInt = sum(isIntRow & isIntCol);
nPyrPyr = sum(~isIntRow & ~isIntCol);
nIntPyr = sum(xor(isIntRow, isIntCol));

counts = [nIntInt, nIntPyr, nPyrPyr];
end
