function runIntPyrLagSkewFDRPermutationTest(saveFile, analysisType, alpha, nNullDraws, flag)
% runs int-vs-pyr-only FDR analysis and skew permutation testing
% for either no-chunk or chunked pairwise saved results

% inputs
%   saveFile      : combined pairwise .mat file
%   analysisType  : "nochunk" or "chunked"
%   alpha         : fdr alpha (default 0.05)
%   nNullDraws    : number of null skew draws (default 100)
%   flag          : passed to david's FDRcutoff (default true)

% outputs
%   saves one .mat file with all results
%   makes plots for each session:
%       1) actual int-vs-pyr significant lag histogram
%       2) pooled null significant lag histogram
%       3) null skew histogram with actual skew marked

% requires david's FDRcutoff.m on the matlab path

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

if exist('FDRcutoff','file') ~= 2
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

    % real int x pyr block
    intRows = 1:nInt;
    pyrCols = nInt + (1:nPyr);

    actual = computeBlockStats( ...
        peakCorr(intRows, pyrCols), ...
        peakLag(intRows, pyrCols), ...
        nullCorr(intRows, pyrCols, :), ...
        alpha, flag);

    fprintf('%s actual: uncorr sig = %d | fdr sig = %d | skew = %.6f\n', ...
        animalID, actual.nSigUncorr, actual.nSigFDR, actual.skew);

    % null skew distribution by random all-vs-all cross-group draws
    nullSkews = nan(nNullDraws,1);
    nullPSkew = nan(nNullDraws,1); %#ok<NASGU>
    nullNSigUncorr = nan(nNullDraws,1);
    nullNSigFDR = nan(nNullDraws,1);
    nullLagCell = cell(nNullDraws,1);

    for r = 1:nNullDraws
        pick = randperm(nAll, nInt + nPyr);   % unique bag of neurons
        grpA = pick(1:nInt);
        grpB = pick(nInt+1:end);

        nullDraw = computeBlockStats( ...
            peakCorr(grpA, grpB), ...
            peakLag(grpA, grpB), ...
            nullCorr(grpA, grpB, :), ...
            alpha, flag);

        nullSkews(r) = nullDraw.skew;
        nullNSigUncorr(r) = nullDraw.nSigUncorr;
        nullNSigFDR(r) = nullDraw.nSigFDR;
        nullLagCell{r} = nullDraw.sigLagVec;
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
    % plot 1: actual significant lag histogram
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

    % plot 2: pooled null significant lag histogram
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

    % plot 3: null skew histogram with actual skew marked
    figure('Name', sprintf('%s skew null distribution', animalID), 'Color', 'w');
    if ~isempty(validNullSkews)
        histogram(validNullSkews, 30, 'FaceAlpha', 0.8, 'EdgeColor', 'none'); hold on;
        xline(actual.skew, 'r-', 'LineWidth', 2);
    else
        histogram(nan, 30, 'FaceAlpha', 0.8, 'EdgeColor', 'none'); hold on;
    end
    xlabel('Skew = (mean - median) / std');
    ylabel('Count');
    title(sprintf('%s skew null distribution | actual=%.4f | p=%.4f', animalID, actual.skew, skewP));
    legend({'Null skews', 'Actual skew'}, 'Location', 'best');
    grid on;

    % save session results
    R = struct();
    R.animalID = animalID;
    R.nInt = nInt;
    R.nPyr = nPyr;
    R.nAll = nAll;

    R.actual = actual;
    R.actual.nPairs = nInt * nPyr;

    R.nullSkews = nullSkews;
    R.validNullSkews = validNullSkews;
    R.nullNSigUncorr = nullNSigUncorr;
    R.nullNSigFDR = nullNSigFDR;
    R.nullLagCell = nullLagCell;
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
function out = computeBlockStats(realBlock, lagBlock, nullBlock, alpha, flag)
% computes pairwise empirical p-values, FDR significance, significant lag skew

[nR, nC, nNull] = size(nullBlock);

realVec = realBlock(:);
lagVec = lagBlock(:);
null2D = reshape(nullBlock, nR*nC, nNull);

pVals = nan(size(realVec));

for i = 1:numel(realVec)
    thisReal = realVec(i);
    thisLag = lagVec(i);
    thisNull = null2D(i,:);
    thisNull = thisNull(~isnan(thisNull) & isfinite(thisNull));

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

sigLagVec = lagVec(sigFDR);
sigLagVec = sigLagVec(~isnan(sigLagVec) & isfinite(sigLagVec));

out = struct();
out.pVals = pVals;
out.fdrCutoff = fdrCut;
out.nValidTests = nnz(validP);
out.nSigUncorr = nnz(sigUncorr);
out.nSigFDR = nnz(sigFDR);
out.sigLagVec = sigLagVec;
out.skew = computeSkew(sigLagVec);

out.pMat = reshape(pVals, nR, nC);
out.sigUncorrMask = reshape(sigUncorr, nR, nC);
out.sigFDRMask = reshape(sigFDR, nR, nC);
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
