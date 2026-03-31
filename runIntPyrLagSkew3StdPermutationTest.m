function runIntPyrLagSkew3StdPermutationTest(saveFile, analysisType, nNullDraws)
% runs int v pyr only 3stdev significance analysis & skew permutation testing for no-chunk & chunked pairwise saved results

% new logic (mimics what i first did way back when):
%   - for each true int x pyr pair, compute null mean and null std from shifted controls
%   - classify pairs significant if real peak correlation > null mean + 3 * null std
%   - compute skew = (mean - median) / std from significant lags

% null test:
%   - sample the same number of unique pairs from all-vs-all upper-triangle
%   - sampled pairs may include int-int, pyr-pyr, and int-pyr
%   - no self-pairs and no duplicate pairs within a draw
%   - compute significant lags and skew using the same 3-std rule
%   - repeat nNullDraws times

% inputs
%   saveFile : combined pairwise .mat file
%   analysisType : "nochunk" or "chunked"
%   nNullDraws : number of null skew draws (default 100)

% output: saves results struct and makes plots per session

if nargin < 3 || isempty(nNullDraws)
    nNullDraws = 100;
end

analysisType = lower(string(analysisType));

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
results.nNullDraws = nNullDraws;
results.sigMethod = 'real_peakCorr_gt_nullMean_plus_3std';
results.sessions = cell(1, nSess);

for s = 1:nSess

    [peakCorr, peakLag, nullCorr, nInt, nPyr, nAll, animalID] = getSessionData(S, s, analysisType);

    fprintf('\n=============================\n');
    fprintf('processing %s (%s)\n', animalID, analysisType);
    fprintf('nInt = %d | nPyr = %d | nAll = %d\n', nInt, nPyr, nAll);
    fprintf('=============================\n');

    % ---------------- actual int x pyr pairs ----------------
    [actualRows, actualCols] = getIntPyrPairs(nInt, nPyr);

    actual = computePairSetStats3Std(peakCorr, peakLag, nullCorr, actualRows, actualCols);

    fprintf('%s actual: nominal pairs = %d | valid tests = %d | sig 3std = %d | skew = %.6f\n', ...
        animalID, numel(actualRows), actual.nValidTests, actual.nSig3Std, actual.skew);

    % ---------------- null draws from all-vs-all upper triangle ----------------
    [poolRows, poolCols] = getAllUpperPairs(nAll);

    nActualNominalPairs = numel(actualRows);

    if nActualNominalPairs > numel(poolRows)
        error('%s: actual nominal pair count exceeds all-vs-all upper-triangle pool size.', animalID);
    end

    nullSkews = nan(nNullDraws,1);
    nullNSig3Std = nan(nNullDraws,1);
    nullNValidTests = nan(nNullDraws,1);
    nullLagCell = cell(nNullDraws,1);
    nullPairTypeCounts = nan(nNullDraws, 3); % [int-int, int-pyr, pyr-pyr]

    for r = 1:nNullDraws
        drawIdx = randperm(numel(poolRows), nActualNominalPairs);
        drawRows = poolRows(drawIdx);
        drawCols = poolCols(drawIdx);

        nullDraw = computePairSetStats3Std(peakCorr, peakLag, nullCorr, drawRows, drawCols);

        nullSkews(r) = nullDraw.skew;
        nullNSig3Std(r) = nullDraw.nSig3Std;
        nullNValidTests(r) = nullDraw.nValidTests;
        nullLagCell{r} = nullDraw.sigLagVec;

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

    figure('Name', sprintf('%s actual int-pyr significant lags (3std)', animalID), 'Color', 'w');
    if ~isempty(actual.sigLagVec)
        lagEdgesActual = makeLagEdges(actual.sigLagVec);
        histogram(actual.sigLagVec, 'BinEdges', lagEdgesActual, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    else
        histogram(nan, 'BinEdges', -0.201:0.001:0.201, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    end
    xlabel('Peak Lag (s)');
    ylabel('Count');
    title(sprintf('%s actual int-pyr 3std significant lags | sig=%d | skew=%.4f', ...
        animalID, actual.nSig3Std, actual.skew));
    grid on;

    pooledNullLags = vertcat(nullLagCell{:});
    pooledNullLags = pooledNullLags(~isnan(pooledNullLags) & isfinite(pooledNullLags));

    figure('Name', sprintf('%s null significant lags (3std)', animalID), 'Color', 'w');
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

    figure('Name', sprintf('%s skew null distribution (3std)', animalID), 'Color', 'w');
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
    R.nullNSig3Std = nullNSig3Std;
    R.nullNValidTests = nullNValidTests;
    R.nullLagCell = nullLagCell;
    R.nullPairTypeCounts = nullPairTypeCounts;
    R.skewP = skewP;

    results.sessions{s} = R;
end

[folderPath, baseName, ~] = fileparts(char(saveFile));
outFile = fullfile(folderPath, sprintf('%s_intPyrSkew3Std_%s.mat', baseName, char(analysisType)));
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
function out = computePairSetStats3Std(realMat, lagMat, nullMat, rows, cols)
% computes pairwise 3-std significance, significant lags, and skew
% for an arbitrary list of selected pairs

nPairs = numel(rows);

realVals = nan(nPairs,1);
lagVals = nan(nPairs,1);
nullMeans = nan(nPairs,1);
nullStds = nan(nPairs,1);
sigThresh = nan(nPairs,1);
sig3Std = false(nPairs,1);

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
        continue;
    end

    thisMean = mean(thisNull, 'omitnan');
    thisStd = std(thisNull, 0, 'omitnan');

    nullMeans(i) = thisMean;
    nullStds(i) = thisStd;
    sigThresh(i) = thisMean + 3*thisStd;

    if isfinite(thisStd) && ~isnan(thisStd) && thisStd > 0
        sig3Std(i) = thisReal > sigThresh(i);
    end
end

validMask = ~isnan(realVals) & isfinite(realVals) & ~isnan(lagVals) & isfinite(lagVals) & ...
            ~isnan(nullMeans) & isfinite(nullMeans) & ~isnan(nullStds) & isfinite(nullStds);

sigLagVec = lagVals(sig3Std);
sigLagVec = sigLagVec(~isnan(sigLagVec) & isfinite(sigLagVec));

out = struct();
out.rows = rows;
out.cols = cols;
out.realVals = realVals;
out.lagVals = lagVals;
out.nullMeans = nullMeans;
out.nullStds = nullStds;
out.sigThresh = sigThresh;
out.nPairsNominal = nPairs;
out.nValidTests = nnz(validMask);
out.nSig3Std = nnz(sig3Std);
out.sigLagVec = sigLagVec;
out.skew = computeSkew(sigLagVec);
out.sig3StdMask = sig3Std;
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
isIntRow = rows <= nInt;
isIntCol = cols <= nInt;

nIntInt = sum(isIntRow & isIntCol);
nPyrPyr = sum(~isIntRow & ~isIntCol);
nIntPyr = sum(xor(isIntRow, isIntCol));

counts = [nIntInt, nIntPyr, nPyrPyr];
end
