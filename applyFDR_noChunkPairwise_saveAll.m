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

S = load(char(saveFile), 'allSessions');
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

    % keep rows where the real value exists
    validReal = ~isnan(realVals) & isfinite(realVals);

    keepLinear = keepLinear(validReal);
    realVals = realVals(validReal);
    nullVals = nullVals(validReal,:);

    pVals = nan(size(realVals));

    for i = 1:numel(realVals)

        thisNull = nullVals(i,:);
        thisNull = thisNull(~isnan(thisNull) & isfinite(thisNull));
        thisReal = realVals(i);

        % require at least a modest number of valid null draws
        if isnan(thisReal) || ~isfinite(thisReal) || numel(thisNull) < 10
            pVals(i) = NaN;
            continue;
        end

        switch tailType

            case "right"
                pVals(i) = (sum(thisNull >= thisReal) + 1) / (numel(thisNull) + 1);

            case "left"
                pVals(i) = (sum(thisNull <= thisReal) + 1) / (numel(thisNull) + 1);

            case "two"

                mu = mean(thisNull, 'omitnan');

                realDev = abs(thisReal - mu);
                nullDev = abs(thisNull - mu);

                pVals(i) = (sum(nullDev >= realDev) + 1) / (numel(thisNull) + 1);

            otherwise
                error('tailType must be "right", "left", or "two".')

        end
    end

    % only pass valid p-values into FDR
    validP = ~isnan(pVals) & isfinite(pVals);

    if any(validP)
        fdrCut = FDRcutoff(pVals(validP), alpha, flag);
    else
        fdrCut = 0;
    end

    sigVec = false(size(pVals));
    sigVec(validP) = pVals(validP) <= fdrCut;

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

    out.nTests = nnz(validP);
    out.nSignificant = nnz(triu(sigMaskFDR,1));
    out.fdrCutoff = fdrCut;

    out.pMat = pMat;
    out.sigMaskFDR = sigMaskFDR;

    out.sigPeakCorrMatAll = sigPeakCorrMatAll;
    out.sigPeakLagMatAll = sigPeakLagMatAll;
    out.sigXcMatAll = sigXcMatAll;

    FDRresults.sessions{s} = out;

    fprintf('%s | valid tests=%d | significant=%d | cutoff=%.6g\n', ...
        animalID, out.nTests, out.nSignificant, out.fdrCutoff);

end

[folderPath, baseName, ~] = fileparts(char(saveFile));
outFile = fullfile(folderPath, [baseName '_FDR.mat']);

save(char(outFile), 'FDRresults', '-v7.3');

fprintf('\nsaved:\n%s\n', outFile);

end
