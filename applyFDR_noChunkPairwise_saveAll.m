function applyFDR_noChunkPairwise_saveAll(saveFile, alpha, flag, tailType)
% applies david's FDRcutoff.m session by session to no-chunk pairwise results
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

S = load(saveFile, 'allSessions');
allSessions = S.allSessions;

nSess = numel(allSessions.sessions);

FDRresults = struct();
FDRresults.sourceFile = saveFile;
FDRresults.alpha = alpha;
FDRresults.flag = flag;
FDRresults.tailType = char(tailType);
FDRresults.sessions = cell(1, nSess);

for s = 1:nSess
    sess = allSessions.sessions(s);

    animalID = regexp(sess.baseDir, 'D\d+', 'match', 'once');
    if isempty(animalID)
        animalID = sprintf('Session%d', s);
    end

    fprintf('\nprocessing %s\n', animalID);

    realMat = sess.peakCorrMatAll;
    lagMat  = sess.peakLagMatAll;
    nullMat = sess.nullCorrMatAllShifts;
    xcMat   = sess.xcMatAll;

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

    sigPeakLagMatAll = lagMat;
    sigPeakLagMatAll(~sigMaskFDR) = NaN;

    sigXcMatAll = xcMat;
    for i = 1:n
        for j = 1:n
            if ~sigMaskFDR(i,j)
                sigXcMatAll(i,j,:) = NaN;
            end
        end
    end

    out = struct();
    out.sessInd = sess.sessInd;
    out.baseDir = sess.baseDir;
    out.combinedFile = sess.combinedFile;
    out.lags = sess.lags;
    out.binSize = sess.binSize;
    out.nInt = sess.nInt;
    out.nPyr = sess.nPyr;
    out.nAll = sess.nAll;
    out.typeVec = sess.typeVec;
    out.numBins = sess.numBins;
    out.xcMatAll = sess.xcMatAll;
    out.peakCorrMatAll = sess.peakCorrMatAll;
    out.peakLagMatAll = sess.peakLagMatAll;
    out.nullCorrMatAllShifts = sess.nullCorrMatAllShifts;
    out.nullCorrPrctile2p5 = sess.nullCorrPrctile2p5;
    out.nullCorrPrctile97p5 = sess.nullCorrPrctile97p5;
    out.nullCorrMean = sess.nullCorrMean;
    out.nullCorrStd = sess.nullCorrStd;
    out.realRowDone = sess.realRowDone;
    out.shiftDone = sess.shiftDone;
    out.realRows = sess.realRows;
    out.shiftNums = sess.shiftNums;
    out.shiftAmtPerNeuronAll = sess.shiftAmtPerNeuronAll;

    out.animalID = animalID;
    out.alpha = alpha;
    out.flag = flag;
    out.tailType = char(tailType);
    out.nTests = numel(pVals);
    out.nSignificant = nnz(triu(sigMaskFDR,1));
    out.fdrCutoff = fdrCut;
    out.pMat = pMat;
    out.sigMaskFDR = sigMaskFDR;
    out.sigPeakCorrMatAll = sigPeakCorrMatAll;
    out.sigPeakLagMatAll = sigPeakLagMatAll;
    out.sigXcMatAll = sigXcMatAll;

    FDRresults.sessions{s} = out;

    fprintf('%s | tests=%d | significant=%d | cutoff=%.6g\n', ...
        animalID, out.nTests, out.nSignificant, out.fdrCutoff);
end

[folderPath, baseName, ~] = fileparts(saveFile);
outFile = fullfile(folderPath, [baseName '_FDR.mat']);
save(outFile, 'FDRresults', '-v7.3');

fprintf('\nsaved:\n%s\n', outFile);
end
