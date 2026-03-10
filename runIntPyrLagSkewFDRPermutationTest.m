function runIntPyrLagSkewFDRPermutationTest(saveFile, analysisType, alpha, nNullDraws, flag)
% runs int-vs-pyr-only FDR analysis and skew permutation testing
% for either no-chunk or chunked pairwise saved results

% actual test:
%   - use only true int x pyr pairs
%   - compute empirical p-values from each pair's shifted null
%   - apply david's FDRcutoff.m
%   - keep FDR-significant lags
%   - compute skew = (mean - median) / std

% null test:
%   - sample the same number of unique pairs from all-vs-all upper-triangle
%   - sampled pairs may include int-int, pyr-pyr, and int-pyr
%   - no self-pairs and no duplicate pairs within a draw
%   - compute p-values, FDR, significant lags, and skew
%   - repeat nNullDraws times

% inputs
%   saveFile      : combined pairwise .mat file
%   analysisType  : "nochunk" or "chunked"
%   alpha         : FDR alpha (default 0.05)
%   nNullDraws    : number of null skew draws (default 100)
%   flag          : passed to david's FDRcutoff (default true)

% output: saves results struct and makes plots per session

if nargin < 3 || isempty(alpha)
    alpha = 0.05;
end
if nargin < 4 || isempty(nNullDraws)
    nNullDraws = 100;
end
if nargin < 5 || isempty(flag)
    flag = true;
end

analysisType = lower(string(analysisType));

if exist('FDRcutoff', 'file') ~= 2
    error('FDRcutoff.m is not on the matlab path.');
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

rng(0); % reproducible null draws

results = struct();
results.sourceFile = saveFile;
results.analysisType = char(analysisType);
results.alpha = alpha;
results.nNullDraws = nNullDraws;
results.flag = flag;
results.sessions = cell(1, nSess);

for s = 1:nSess

    [peakCorr, peakLag, nullCorr, nInt, nPyr, nAll, animalID] = getSessionData(S, s, analysisType);

    fprintf('\n=============================\n');
    fprintf('processing %s (%s)\n', animalID, analysisType);
    fprintf('nInt = %d | nPyr = %d | nAll = %d\n', nInt, nPyr, nAll);
    fprintf('=============================\n');

    % ---------------- actual int x pyr pairs ----------------
    [actualRows, actualCols] = getIntPyrPairs(nInt, nPyr);

    actual = computePairSetStats(peakCorr, peakLag, nullCorr, actualRows, actualCols, alpha, flag);

    fprintf('%s actual: nominal pairs = %d | valid tests = %d | uncorr sig = %d | fdr sig = %d | skew = %.6f\n', ...
        animalID, numel(actualRows), actual.nValidTests, actual.nSigUncorr, actual.nSigFDR, actual.skew);

    % ---------------- null draws from all-vs-all upper triangle ----------------
    [poolRows, poolCols] = getAllUpperPairs(nAll);

    % keep only unique upper-triangle pairs by construction
    nActualNominalPairs = numel(actualRows);

    if nActualNominalPairs > numel(poolRows)
        error('%s: actual nominal pair count exceeds all-vs-all upper-triangle pool size.', animalID);
    end

    nullSkews = nan(nNullDraws,1);
    nullNSigUncorr = nan(nNullDraws,1);
    nullNSigFDR = nan(nNullDraws,1);
    nullNValidTests = nan(nNullDraws,1);
    nullLagCell = cell(nNullDraws,1);
    nullFdrCutoff = nan(nNullDraws,1);
    nullPairTypeCounts = nan(nNullDraws, 3); % [int-int, int-pyr, pyr-pyr]

    for r = 1:nNullDraws
        drawIdx = randperm(numel(poolRows), nActualNominalPairs);  % no duplicates, no self-pairs
        drawRows = poolRows(drawIdx);
        drawCols = poolCols(drawIdx);

        nullDraw = computePairSetStats(peakCorr, peakLag, nullCorr, drawRows, drawCols, alpha, flag);

        nullSkews(r) = nullDraw.skew;
        nullNSigUncorr(r) = nullDraw.nSigUncorr;
        nullNSigFDR(r) = nullDraw.nSigFDR;
        nullNValidTests(r) = nullDraw.nValidTests;
        nullLagCell{r} = nullDraw.sigLagVec;
        nullFdrCutoff(r) = nullDraw.fdrCutoff;

        nullPairTypeCounts(r,:) = classifySampledPairs(drawRows, drawCols, nInt);
    end

    validNullSkews = nullSkews(~isnan(nullSkews) & isfinite(nullSkews));

    if isempty(validNullSkews) || isnan(actual.skew) || ~isfinite(actual.skew)
        skewP = NaN;
    else
        skewP = (sum(validNullSkews >= actual.skew) + 1) / (numel(validNullSkews) + 1);
    end

    fprintf('%s skew test: actual skew = %.6f | null n = %d | skew p = %.6f\n', ...
        animalID, actual.skew, numel(validNullSkews), skewP);

    % ---------------- plots ----------------

    % actual significant lag histogram
    figure('Name', sprintf('%s actual int-pyr significant lags', animalID), 'Color', 'w');
    if ~isempty(actual.sigLagVec)
        lagEdgesActual = makeLagEdges(actual.sigLagVec);
        histogram(actual.sigLagVec, 'BinEdges', lagEdgesActual, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    else
        histogram(nan, 'BinEdges', -0.201:0.001:0.201, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    end
    xlabel('Peak Lag (s)');
    ylabel('Count');
    title(sprintf('%s actual int-pyr FDR-significant lags | uncorr=%d | fdr=%d | skew=%.4f', ...
        animalID, actual.nSigUncorr, actual.nSigFDR, actual.skew));
    grid on;

    % pooled null significant lag histogram
    pooledNullLags = vertcat(nullLagCell{:});
    pooledNullLags = pooledNullLags(~isnan(pooledNullLags) & isfinite(pooledNullLags));

    figure('Name', sprintf('%s null significant lags', animalID), 'Color', 'w');
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

    % null skew histogram with actual skew marked
    figure('Name', sprintf('%s skew null distribution', animalID), 'Color', 'w');
    if ~isempty(validNullSkews)
        histogram(validNullSkews, 30, 'FaceAlpha', 0.8, 'EdgeColor', 'none'); hold on;
        xline(actual.skew, 'r-', 'LineWidth', 2);
        legend({'Null skews', 'Actual skew'}, 'Location', 'best');
    else
        histogram(nan, 30, 'FaceAlpha', 0.8, 'EdgeColor', 'none'); hold on;
        legend({'Null skews'}, 'Location', 'best');
    end
    xlabel('Skew = (mean - median) / std');
    ylabel('Count');
    title(sprintf('%s skew null distribution | actual=%.4f | p=%.4f', animalID, actual.skew, skewP));
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

    R.nullSkews = nullSkews;
    R.validNullSkews = validNullSkews;
    R.nullNSigUncorr = nullNSigUncorr;
    R.nullNSigFDR = nullNSigFDR;
    R.nullNValidTests = nullNValidTests;
    R.nullLagCell = nullLagCell;
    R.nullFdrCutoff = nullFdrCutoff;
    R.nullPairTypeCounts = nullPairTypeCounts;
    R.skewP = skewP;

    results.sessions{s} = R;
end

[folderPath, baseName, ~] = fileparts(char(saveFile));
outFile = fullfile(folderPath, sprintf('%s_intPyrSkewFDR_%s.mat', baseName, char(analysisType)));
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
% true int x pyr block: rows are int, cols are pyr
[rr, cc] = ndgrid(1:nInt, nInt + (1:nPyr));
rows = rr(:);
cols = cc(:);
end

% ============================================================
function [rows, cols] = getAllUpperPairs(nAll)
% all-vs-all unique upper-triangle pairs, excluding self-pairs
mask = triu(true(nAll), 1);
[rows, cols] = find(mask);
end

% ============================================================
function out = computePairSetStats(realMat, lagMat, nullMat, rows, cols, alpha, flag)
% computes pairwise empirical p-values, FDR significance, and skew
% for an arbitrary list of selected pairs

nPairs = numel(rows);

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
    lagVals(i) = thisLag;

    if isnan(thisReal) || ~isfinite(thisReal) || isnan(thisLag) || ~isfinite(thisLag) || numel(thisNull) < 10
        pVals(i) = NaN;
        continue;
    end

    % right-tailed empirical p-value
    pVals(i) = (sum(thisNull >= thisReal) + 1) / (numel(thisNull) + 1);
end

validP = ~isnan(pVals) & isfinite(pVals);

sigUncorr = false(size(pVals));
sigUncorr(validP) = pVals(validP) <= alpha;

if any(validP)
    fdrCut = FDRcutoff(pVals(validP), alpha, flag);
else
    fdrCut = 0;
end

sigFDR = false(size(pVals));
sigFDR(validP) = pVals(validP) <= fdrCut;

sigLagVec = lagVals(sigFDR);
sigLagVec = sigLagVec(~isnan(sigLagVec) & isfinite(sigLagVec));

out = struct();
out.rows = rows;
out.cols = cols;
out.realVals = realVals;
out.lagVals = lagVals;
out.pVals = pVals;
out.fdrCutoff = fdrCut;
out.nPairsNominal = nPairs;
out.nValidTests = nnz(validP);
out.nSigUncorr = nnz(sigUncorr);
out.nSigFDR = nnz(sigFDR);
out.sigLagVec = sigLagVec;
out.skew = computeSkew(sigLagVec);
out.sigUncorrMask = sigUncorr;
out.sigFDRMask = sigFDR;
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
% returns [int-int, int-pyr, pyr-pyr] counts for one sampled null draw

isIntRow = rows <= nInt;
isIntCol = cols <= nInt;

nIntInt = sum(isIntRow & isIntCol);
nPyrPyr = sum(~isIntRow & ~isIntCol);
nIntPyr = sum(xor(isIntRow, isIntCol));

counts = [nIntInt, nIntPyr, nPyrPyr];
end
