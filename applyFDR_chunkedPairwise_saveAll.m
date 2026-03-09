function applyFDR_chunkedPairwise_saveAll(saveFile, alpha, flag, tailType)
% applies david's FDRcutoff.m session by session to chunked pairwise results
% and saves all original variables again plus significant-only versions

if nargin < 2 || isempty(alpha)
    alpha = 0.05;
end
if nargin < 3 || isempty(flag)
    flag = true;
end
if nargin < 4 || isempty(tailType)
    tailType = "right";
end
tailType = lower(string(tailType));

if exist('FDRcutoff','file') ~= 2
    error('FDRcutoff.m is not on the matlab path.');
end

S = load(saveFile);

nSess = numel(S.baseDirs);

FDRresults = struct();
FDRresults.sourceFile = saveFile;
FDRresults.alpha = alpha;
FDRresults.flag = flag;
FDRresults.tailType = char(tailType);
FDRresults.baseDirs = S.baseDirs;
FDRresults.sessNames = S.sessNames;
FDRresults.sessions = repmat(struct(), 1, nSess);

for s = 1:nSess
    baseDir = S.baseDirs{s};

    animalID = regexp(baseDir, 'D\d+', 'match', 'once');
    if isempty(animalID)
        animalID = sprintf('Session%d', s);
    end

    fprintf('\nprocessing %s\n', animalID);

    realMat = S.all_peakCorrMat_all{s};
    lagMat  = S.all_peakLagSecMat_all{s};
    nullMat = S.all_nullCorrMat_allShifts{s};
    xcMat   = S.all_xcMat_all{s};

    n = size(realMat,1);
    upperMask = triu(true(n),1);
    keepLinear = find(upperMask);

    realVals = realMat(keepLinear);

    nNull = size(nullMat,3);
    nullVals = nan(numel(keepLinear), nNull);
    for k = 1:nNull
        tmp = nullMat(:,:,k);
        nullVals(:,k) = tmp(keepLinear);
    end

    valid = ~isnan(realVals) & all(~isnan(nullVals),2);
    keepLinear = keepLinear(valid);
    realVals = realVals(valid);
    nullVals = nullVals(valid,:);

    pVals = nan(size(realVals));

    for i = 1:numel(realVals)
        switch tailType
            case "right"
                pVals(i) = (sum(nullVals(i,:) >= realVals(i)) + 1) / (nNull + 1);
            case "left"
                pVals(i) = (sum(nullVals(i,:) <= realVals(i)) + 1) / (nNull + 1);
            case "two"
                mu = mean(nullVals(i,:), 'omitnan');
                realDev = abs(realVals(i) - mu);
                nullDev = abs(nullVals(i,:) - mu);
                pVals(i) = (sum(nullDev >= realDev) + 1) / (nNull + 1);
            otherwise
                error('tailType must be "right", "left", or "two".');
        end
    end

    fdrCut = FDRcutoff(pVals, alpha, flag);
    sigVec = pVals <= fdrCut;

    pMat = nan(n);
    sigMaskFDR = false(n);

    pMat(keepLinear) = pVals;
    sigMaskFDR(keepLinear) = sigVec;

    pMat = pMat + pMat.';
    sigMaskFDR = sigMaskFDR | sigMaskFDR.';

    sigPeakCorrMatAll = realMat;
    sigPeakCorrMatAll(~sigMaskFDR) = NaN;

    sigPeakLagSecMat_all = lagMat;
    sigPeakLagSecMat_all(~sigMaskFDR) = NaN;

    sigXcMatAll = xcMat;
    for i = 1:n
        for j = 1:n
            if ~sigMaskFDR(i,j)
                sigXcMatAll(i,j,:) = NaN;
            end
        end
    end

    out = struct();
    out.baseDir = baseDir;
    out.animalID = animalID;

    % save originals again
    out.binSize = S.all_binSize{s};
    out.channelsToUse = S.all_channelsToUse{s};
    out.chunkHalf = S.all_chunkHalf{s};
    out.doBaselineNorm = S.all_doBaselineNorm{s};
    out.intIdx = S.all_intIdx{s};
    out.lags = S.all_lags{s};
    out.nAll = S.all_nAll(s);
    out.nInt_ref = S.all_nInt_ref(s);
    out.nPyr_ref = S.all_nPyr_ref(s);
    out.neuronType = S.all_neuronType{s};
    out.nullCorrMat_allShifts = nullMat;
    out.peakCorrMat_all = realMat;
    out.peakLagSecMat_all = lagMat;
    out.pyrIdx = S.all_pyrIdx{s};
    out.runtimeRealSec = S.all_runtimeRealSec(s);
    out.runtimeShiftSec = S.all_runtimeShiftSec{s};
    out.tAxis = S.all_tAxis{s};
    out.xcMat_all = xcMat;
    out.realFound = S.realFound(s);
    out.shiftFound = S.shiftFound(s,:);
    out.sessName = S.sessNames{s};

    % save FDR results
    out.alpha = alpha;
    out.flag = flag;
    out.tailType = char(tailType);
    out.nTests = numel(pVals);
    out.nSignificant = nnz(triu(sigMaskFDR,1));
    out.fdrCutoff = fdrCut;
    out.pMat = pMat;
    out.sigMaskFDR = sigMaskFDR;
    out.sigPeakCorrMat_all = sigPeakCorrMatAll;
    out.sigPeakLagSecMat_all = sigPeakLagSecMat_all;
    out.sigXcMat_all = sigXcMatAll;

    FDRresults.sessions(s) = out;

    fprintf('%s | tests=%d | significant=%d | cutoff=%.6g\n', ...
        animalID, out.nTests, out.nSignificant, out.fdrCutoff);
end

[folderPath, baseName, ~] = fileparts(saveFile);
outFile = fullfile(folderPath, [baseName '_FDR.mat']);
save(outFile, 'FDRresults', '-v7.3');

fprintf('\nsaved:\n%s\n', outFile);
end
