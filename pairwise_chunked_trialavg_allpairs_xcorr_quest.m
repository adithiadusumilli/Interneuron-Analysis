% pairwise_chunked_trialavg_allpairs_xcorr_quest.m
% jobInd = 0 -> real all-pairs trial-averaged xcorr
% jobInd = 1:100 -> random ordered all-vs-all null peak lag draws

clc;

baseDirs = {
    '/home/asa7288/Transfer/D026', ...
    '/home/asa7288/Transfer/D020', ...
    '/home/asa7288/Transfer/D024', ...
    '/home/asa7288/Transfer/D043'
};

sessInd = 4; 
baseDirs = baseDirs(sessInd);

channelsToUse  = 1:4;
chunkHalf = 300;
maxLagSecs = 0.2;
doBaselineNorm = true;

assert(exist('jobInd','var')==1, 'define jobInd before running');
assert(jobInd >= 0 && jobInd <= 100, 'jobInd must be 0:100');

rng(jobInd);

for iDir = 1:numel(baseDirs)

    baseDir = baseDirs{iDir};
    fprintf('\nprocessing %s | sessInd=%d | jobInd=%d\n', baseDir, sessInd, jobInd);

    S = load(fullfile(baseDir, 'EMG_Neural_AllChannels.mat'), 'pyrCxWinCell', 'intCxWinCell', 'tAxis');

    tAxis = S.tAxis(:)';
    binSize = 0.001;
    T = numel(tAxis);

    [~, cIdx] = min(abs(tAxis - 0));
    chunkStart = cIdx - chunkHalf;
    chunkEnd = cIdx + chunkHalf;

    reqMaxLagBins = round(maxLagSecs / binSize);
    maxLagLeft = chunkStart - 1;
    maxLagRight = T - chunkEnd;
    maxLagBins = min([reqMaxLagBins, maxLagLeft, maxLagRight]);
    lags = -maxLagBins:maxLagBins;

    pyrAll = [];
    intAll = [];
    nPyr_ref = nan;
    nInt_ref = nan;

    for ci = 1:numel(channelsToUse)
        ch = channelsToUse(ci);

        pyrWin = S.pyrCxWinCell{ch};
        intWin = S.intCxWinCell{ch};

        if isempty(pyrWin) || isempty(intWin)
            continue;
        end

        nEvt = min(size(pyrWin,1), size(intWin,1));
        pyrWin = pyrWin(1:nEvt,:,:);
        intWin = intWin(1:nEvt,:,:);

        if isnan(nPyr_ref), nPyr_ref = size(pyrWin,2); end
        if isnan(nInt_ref), nInt_ref = size(intWin,2); end

        nPyr_ref = min(nPyr_ref, size(pyrWin,2));
        nInt_ref = min(nInt_ref, size(intWin,2));

        pyrWin = pyrWin(:,1:nPyr_ref,:);
        intWin = intWin(:,1:nInt_ref,:);

        if doBaselineNorm
            pyrWin = subtractTrialBaseline3d(pyrWin, tAxis, -500, -450);
            intWin = subtractTrialBaseline3d(intWin, tAxis, -500, -450);
        end

        pyrAll = cat(1, pyrAll, pyrWin);
        intAll = cat(1, intAll, intWin);
    end

    if isempty(pyrAll) || isempty(intAll)
        warning('no valid data for %s', baseDir);
        continue;
    end

    nEvt = min(size(pyrAll,1), size(intAll,1));
    pyrAll = pyrAll(1:nEvt,:,:);
    intAll = intAll(1:nEvt,:,:);

    pyrAvg = squeeze(mean(pyrAll, 1, 'omitnan'));
    intAvg = squeeze(mean(intAll, 1, 'omitnan'));

    allAvg = cat(1, pyrAvg, intAvg);
    nAll = size(allAvg,1);

    neuronType = [zeros(1,nPyr_ref), ones(1,nInt_ref)];
    pyrIdx = 1:nPyr_ref;
    intIdx = (nPyr_ref+1):nAll;

    outDir = fullfile(baseDir, 'quest_runs');
    if ~exist(outDir,'dir'), mkdir(outDir); end

    tic;

    if jobInd == 0

        [xcMat_all, peakCorrMat_all, peakLagSecMat_all] = lagSweepPairwise_fullxc_trialavg_allpairs(allAvg, lags, chunkStart, chunkEnd, T, binSize);

        runtimeSec = toc;

        outFile = fullfile(outDir, ...
            sprintf('pairwiseChunkedTrialAvgXCorr_ALLPAIRS_unshifted_sess%02d_fullxc.mat', sessInd));

        save(outFile, ...
            'lags','binSize','tAxis', ...
            'xcMat_all','peakCorrMat_all','peakLagSecMat_all', ...
            'chunkHalf','channelsToUse','doBaselineNorm', ...
            'jobInd','baseDir','sessInd', ...
            'nPyr_ref','nInt_ref','nAll','neuronType','pyrIdx','intIdx', ...
            'runtimeSec','-v7.3');

    else

        nDrawPairs = nAll * (nAll - 1);

        [nullPeakLagVec, nullPeakCorrVec, nullPairI, nullPairJ, nullPairType] = ...
            randomOrderedAllPairsNull_trialavg( ...
                allAvg, neuronType, nDrawPairs, lags, chunkStart, chunkEnd, T, binSize);

        runtimeSec = toc;

        outFile = fullfile(outDir, ...
            sprintf('pairwiseChunkedTrialAvgXCorr_ALLPAIRS_sess%02d_permNull_%03d.mat', sessInd, jobInd));

        save(outFile, ...
            'nullPeakLagVec','nullPeakCorrVec','nullPairI','nullPairJ','nullPairType', ...
            'lags','binSize','tAxis', ...
            'chunkHalf','channelsToUse','doBaselineNorm', ...
            'jobInd','baseDir','sessInd', ...
            'nPyr_ref','nInt_ref','nAll','neuronType','pyrIdx','intIdx', ...
            'runtimeSec','-v7.3');
    end

    fprintf('saved: %s\n', outFile);
end

function Wout = subtractTrialBaseline3d(Win, tAxis, tStart, tEnd)
    [~, iStart] = min(abs(tAxis - tStart));
    [~, iEnd]   = min(abs(tAxis - tEnd));
    if iEnd < iStart
        tmp = iStart; iStart = iEnd; iEnd = tmp;
    end
    base = mean(Win(:,:,iStart:iEnd), 3, 'omitnan');
    Wout = Win - base;
end

function [xcMat, peakCorrMat, peakLagSecMat] = lagSweepPairwise_fullxc_trialavg_allpairs(allAvg, lags, chunkStart, chunkEnd, T, binSize)

    nAll = size(allAvg,1);
    nL = numel(lags);

    xcMat = nan(nAll, nAll, nL);
    peakCorrMat = nan(nAll, nAll);
    peakLagSecMat = nan(nAll, nAll);

    for i = 1:nAll
        for j = (i+1):nAll

            [xc, peakLagSec, peakCorr] = onePairLagSweep( ...
                allAvg(i,:), allAvg(j,:), lags, chunkStart, chunkEnd, T, binSize);

            xcMat(i,j,:) = xc;
            peakCorrMat(i,j) = peakCorr;
            peakLagSecMat(i,j) = peakLagSec;
        end
    end
end

function [nullPeakLagVec, nullPeakCorrVec, nullPairI, nullPairJ, nullPairType] = randomOrderedAllPairsNull_trialavg(allAvg, neuronType, nDrawPairs, lags, chunkStart, chunkEnd, T, binSize)

    nAll = size(allAvg,1);

    nullPeakLagVec = nan(nDrawPairs,1);
    nullPeakCorrVec = nan(nDrawPairs,1);
    nullPairI = nan(nDrawPairs,1);
    nullPairJ = nan(nDrawPairs,1);
    nullPairType = strings(nDrawPairs,1);

    for d = 1:nDrawPairs

        i = randi(nAll);
        j = randi(nAll);

        while j == i
            j = randi(nAll);
        end

        [~, peakLagSec, peakCorr] = onePairLagSweep( ...
            allAvg(i,:), allAvg(j,:), lags, chunkStart, chunkEnd, T, binSize);

        nullPeakLagVec(d) = peakLagSec;
        nullPeakCorrVec(d) = peakCorr;
        nullPairI(d) = i;
        nullPairJ(d) = j;

        if neuronType(i)==0 && neuronType(j)==0
            nullPairType(d) = "pyr-pyr";
        elseif neuronType(i)==1 && neuronType(j)==1
            nullPairType(d) = "int-int";
        elseif neuronType(i)==0 && neuronType(j)==1
            nullPairType(d) = "pyr-int";
        elseif neuronType(i)==1 && neuronType(j)==0
            nullPairType(d) = "int-pyr";
        end
    end
end

function [xc, peakLagSec, peakCorr] = onePairLagSweep(aTrace, bTrace, lags, chunkStart, chunkEnd, T, binSize)

    xc = nan(1,numel(lags));

    for iL = 1:numel(lags)
        L = lags(iL);

        bStart = chunkStart - L;
        bEnd   = chunkEnd   - L;

        if bStart < 1 || bEnd > T
            continue;
        end

        aseg = aTrace(chunkStart:chunkEnd);
        bseg = bTrace(bStart:bEnd);

        valid = ~isnan(aseg) & ~isnan(bseg);

        if nnz(valid) > 10
            xc(iL) = corr(aseg(valid)', bseg(valid)');
        end
    end

    [peakCorr, peakIdx] = max(xc);
    peakLagSec = lags(peakIdx) * binSize;
end
