% chunked trial-averaged population cross-correlation
% averages across neurons, baseline-subtracts each event/window,
% averages across trials/windows, then cross-correlates population averages
%
% jobInd meaning:
%   0         -> unpermuted real
%   1:100     -> permutation jobs
%   101:200   -> shifted-control jobs (shift #1..100)
%
% permutation update:
%   for perm jobs, rerun the permutation until the permuted peak correlation
%   is >= the 95th percentile of shifted-control zero-lag correlations
%   for that animal/session, or until maxPermTries is reached

clc;

%% ---------------- settings ----------------
baseDirs = {
    '/home/asa7288/Transfer/D026', ...
    '/home/asa7288/Transfer/D020', ...
    '/home/asa7288/Transfer/D024', ...
    '/home/asa7288/Transfer/D043'
};

rng(jobInd)

channelsToUse  = 1:4;
chunkHalf      = 300;   % -300 to +300 ms
maxLagSecs     = 0.2;   % -200 to +200 ms
doBaselineNorm = true;

maxPermTries = 1000;

%% ---------------- job type ----------------
if jobInd == 0
    jobType = "real";
    permInd = 0;
    shiftInd = 0;
elseif jobInd >= 1 && jobInd <= 100
    jobType = "perm";
    permInd = jobInd;
    shiftInd = 0;
elseif jobInd >= 101 && jobInd <= 200
    jobType = "shift";
    permInd = 0;
    shiftInd = jobInd - 100;
else
    error('jobInd must be 0, 1:100, or 101:200');
end

for iDir = 1:numel(baseDirs)
    baseDir = baseDirs{iDir};
    fprintf('\nprocessing cortex session %d: %s (job %d, type=%s)\n', iDir, baseDir, jobInd, jobType);

    if jobType == "shift"
        S = load(fullfile(baseDir, 'EMG_Neural_AllChannels.mat'), ...
            'pyrCxWinCell','intCxWinCell','tAxis','validTransitionsNeurShiftedCell');
    else
        S = load(fullfile(baseDir, 'EMG_Neural_AllChannels.mat'), ...
            'pyrCxWinCell','intCxWinCell','tAxis');
    end

    tAxis = S.tAxis(:)';
    binSize = 0.001;
    T = numel(tAxis);

    [~, cIdx] = min(abs(tAxis - 0));

    pyrStart = cIdx - chunkHalf;
    pyrEnd   = cIdx + chunkHalf;

    reqMaxLagBins = round(maxLagSecs / binSize);
    maxLagLeft  = pyrStart - 1;
    maxLagRight = T - pyrEnd;
    maxLagBins  = min([reqMaxLagBins, maxLagLeft, maxLagRight]);
    lags = -maxLagBins:maxLagBins;

    %% ----------- load shift-control upper bound for permutation jobs -----------
    shiftCorrUpper95 = NaN;

    if jobType == "perm"
        shiftCorrUpper95 = getShiftCorrUpper95(baseDir);

        if isnan(shiftCorrUpper95)
            warning('could not compute shiftCorrUpper95 for %s; permutation will run once without rerun filtering.', baseDir);
        else
            fprintf('shiftCorrUpper95 for %s = %.6f\n', baseDir, shiftCorrUpper95);
        end
    end

    %% ----------- compute metric -----------
    tic

    if jobType == "perm"
        tryCount = 0;
        acceptedPerm = false;

        while ~acceptedPerm && tryCount < maxPermTries
            tryCount = tryCount + 1;

            [pyr_byCh, int_byCh] = buildEventsByChannel( ...
                S, channelsToUse, tAxis, doBaselineNorm, jobType);

            [pyrAvgTrace, intAvgTrace, nTrialsUsed] = averageAcrossTrialsAndChannels(pyr_byCh, int_byCh);

            if isempty(pyrAvgTrace) || isempty(intAvgTrace)
                warning('no usable averaged traces for %s, skipping', baseDir);
                break;
            end

            [xc, peakLagSec] = lagSweepTrialAverage(pyrAvgTrace, intAvgTrace, lags, pyrStart, pyrEnd, T, binSize);
            permPeakCorr = max(xc, [], 'omitnan');

            if isnan(shiftCorrUpper95) || permPeakCorr >= shiftCorrUpper95
                acceptedPerm = true;
            end
        end

        runtimeSec = toc;

        fprintf('perm done: accepted=%d | tries=%d | permPeakCorr=%.6f | shiftCorrUpper95=%.6f | time=%.3f s\n', ...
            acceptedPerm, tryCount, permPeakCorr, shiftCorrUpper95, runtimeSec);

    elseif jobType == "shift"
        [pyr_byCh, int_byCh] = buildShiftEventsByChannel( ...
            S, baseDir, channelsToUse, tAxis, doBaselineNorm, shiftInd);

        [pyrAvgTrace, intAvgTrace, nTrialsUsed] = averageAcrossTrialsAndChannels(pyr_byCh, int_byCh);

        if isempty(pyrAvgTrace) || isempty(intAvgTrace)
            warning('no usable averaged traces for %s, skipping', baseDir);
            continue;
        end

        xcZeroLag = zeroLagTrialAverage(pyrAvgTrace, intAvgTrace, pyrStart, pyrEnd);
        runtimeSec = toc;
        disp(['shifted trial-avg control done, time: ' num2str(runtimeSec) ' s'])

        tryCount = NaN;
        acceptedPerm = NaN;
        permPeakCorr = NaN;

    else
        [pyr_byCh, int_byCh] = buildEventsByChannel( ...
            S, channelsToUse, tAxis, doBaselineNorm, jobType);

        [pyrAvgTrace, intAvgTrace, nTrialsUsed] = averageAcrossTrialsAndChannels(pyr_byCh, int_byCh);

        if isempty(pyrAvgTrace) || isempty(intAvgTrace)
            warning('no usable averaged traces for %s, skipping', baseDir);
            continue;
        end

        [xc, peakLagSec] = lagSweepTrialAverage(pyrAvgTrace, intAvgTrace, lags, pyrStart, pyrEnd, T, binSize);
        runtimeSec = toc;
        disp(['trial-averaged cross-correlation analysis done, time: ' num2str(runtimeSec) ' s'])

        tryCount = NaN;
        acceptedPerm = NaN;
        permPeakCorr = NaN;
    end

    %% ----------- save results -----------
    outDir = fullfile(baseDir, 'quest_runs');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    if jobType == "real"
        outFile = fullfile(outDir, 'concatCrossCorr_trialavg_chunked_popavg_unperm.mat');
        save(outFile, ...
            'lags','binSize','xc','peakLagSec','chunkHalf','channelsToUse', ...
            'doBaselineNorm','jobInd','permInd','shiftInd','jobType','baseDir', ...
            'runtimeSec','pyrAvgTrace','intAvgTrace','nTrialsUsed', ...
            'tryCount','acceptedPerm','permPeakCorr','shiftCorrUpper95');
    elseif jobType == "perm"
        outFile = fullfile(outDir, sprintf('concatCrossCorr_trialavg_chunked_popavg_perm_%03d.mat', permInd));
        save(outFile, ...
            'lags','binSize','xc','peakLagSec','chunkHalf','channelsToUse', ...
            'doBaselineNorm','jobInd','permInd','shiftInd','jobType','baseDir', ...
            'runtimeSec','pyrAvgTrace','intAvgTrace','nTrialsUsed', ...
            'tryCount','acceptedPerm','permPeakCorr','shiftCorrUpper95','maxPermTries');
    else
        outFile = fullfile(outDir, sprintf('concatCrossCorr_trialavg_chunked_popavg_shift_%03d_zerolag.mat', shiftInd));
        save(outFile, ...
            'binSize','xcZeroLag','chunkHalf','channelsToUse', ...
            'doBaselineNorm','jobInd','permInd','shiftInd','jobType','baseDir', ...
            'runtimeSec','pyrAvgTrace','intAvgTrace','nTrialsUsed');
    end

    fprintf('saved: %s\n', outFile);
end

%% ===================== helpers =====================

function [pyr_byCh, int_byCh] = buildEventsByChannel(S, channelsToUse, tAxis, doBaselineNorm, jobType)
    pyr_byCh = cell(1, numel(channelsToUse));
    int_byCh = cell(1, numel(channelsToUse));

    for ci = 1:numel(channelsToUse)
        ch = channelsToUse(ci);

        pyrWin = S.pyrCxWinCell{ch};
        intWin = S.intCxWinCell{ch};

        if isempty(pyrWin) || isempty(intWin)
            continue;
        end

        if jobType == "perm"
            [nEvtP, nPyr, ~] = size(pyrWin);
            [nEvtI, nInt, ~] = size(intWin);
            nEvt = min(nEvtP, nEvtI);

            if nPyr + nInt < 2 || nEvt < 1
                continue;
            end

            pooledWin = cat(2, pyrWin(1:nEvt,:,:), intWin(1:nEvt,:,:));
            totN = nPyr + nInt;

            idx = randperm(totN);
            idxP = idx(1:nPyr);
            idxI = idx(nPyr+1:end);

            pyrEvt = squeeze(mean(pooledWin(:, idxP, :), 2, 'omitnan'));
            intEvt = squeeze(mean(pooledWin(:, idxI, :), 2, 'omitnan'));
        else
            pyrEvt = squeeze(mean(pyrWin, 2, 'omitnan'));
            intEvt = squeeze(mean(intWin, 2, 'omitnan'));
        end

        if doBaselineNorm
            pyrEvt = subtractTrialBaseline(pyrEvt, tAxis, -500, -450);
            intEvt = subtractTrialBaseline(intEvt, tAxis, -500, -450);
        end

        pyr_byCh{ci} = pyrEvt;
        int_byCh{ci} = intEvt;
    end
end

function [pyr_byCh, int_byCh] = buildShiftEventsByChannel(S, baseDir, channelsToUse, tAxis, doBaselineNorm, shiftInd)
    pyr_byCh = cell(1, numel(channelsToUse));
    int_byCh = cell(1, numel(channelsToUse));

    [cortexIntFR, okI] = loadCortexIntFR(baseDir);

    if ~okI
        warning('failed to load cortex interneuron FR for %s', baseDir);
        return;
    end

    shiftCol = shiftInd + 1;

    for ci = 1:numel(channelsToUse)
        ch = channelsToUse(ci);

        pyrWin = S.pyrCxWinCell{ch};

        if isempty(pyrWin)
            continue;
        end

        pyrEvt = squeeze(mean(pyrWin, 2, 'omitnan'));

        idx1k = S.validTransitionsNeurShiftedCell{ch, shiftCol};
        if isempty(idx1k)
            continue;
        end

        intWinShift = extractWindowsFromFR(cortexIntFR, idx1k, tAxis);
        intEvt = squeeze(mean(intWinShift, 2, 'omitnan'));

        nEvt = min(size(pyrEvt,1), size(intEvt,1));
        pyrEvt = pyrEvt(1:nEvt,:);
        intEvt = intEvt(1:nEvt,:);

        if doBaselineNorm
            pyrEvt = subtractTrialBaseline(pyrEvt, tAxis, -500, -450);
            intEvt = subtractTrialBaseline(intEvt, tAxis, -500, -450);
        end

        pyr_byCh{ci} = pyrEvt;
        int_byCh{ci} = intEvt;
    end
end

function shiftCorrUpper95 = getShiftCorrUpper95(baseDir)
    outDir = fullfile(baseDir, 'quest_runs');
    shiftFiles = dir(fullfile(outDir, 'concatCrossCorr_trialavg_chunked_popavg_shift_*_zerolag.mat'));

    if isempty(shiftFiles)
        shiftCorrUpper95 = NaN;
        return;
    end

    shiftVals = nan(numel(shiftFiles),1);

    for k = 1:numel(shiftFiles)
        D = load(fullfile(outDir, shiftFiles(k).name), 'xcZeroLag');

        if isfield(D, 'xcZeroLag')
            shiftVals(k) = D.xcZeroLag;
        end
    end

    shiftVals = shiftVals(~isnan(shiftVals) & isfinite(shiftVals));

    if isempty(shiftVals)
        shiftCorrUpper95 = NaN;
    else
        shiftCorrUpper95 = prctile(shiftVals, 95);
    end
end

function Eout = subtractTrialBaseline(Ein, tAxis, tStart, tEnd)
    [~, iStart] = min(abs(tAxis - tStart));
    [~, iEnd]   = min(abs(tAxis - tEnd));

    if iEnd < iStart
        tmp = iStart;
        iStart = iEnd;
        iEnd = tmp;
    end

    baselines = mean(Ein(:, iStart:iEnd), 2, 'omitnan');
    Eout = Ein - baselines;
end

function [pyrAvgTrace, intAvgTrace, nTrialsUsed] = averageAcrossTrialsAndChannels(pyr_byCh, int_byCh)
    pyrAll = [];
    intAll = [];

    for ci = 1:numel(pyr_byCh)
        pEvt = pyr_byCh{ci};
        iEvt = int_byCh{ci};

        if isempty(pEvt) || isempty(iEvt)
            continue;
        end

        nEvt = min(size(pEvt,1), size(iEvt,1));
        if nEvt < 1
            continue;
        end

        pyrAll = cat(1, pyrAll, pEvt(1:nEvt,:));
        intAll = cat(1, intAll, iEvt(1:nEvt,:));
    end

    nTrialsUsed = size(pyrAll, 1);

    if nTrialsUsed < 1
        pyrAvgTrace = [];
        intAvgTrace = [];
        return;
    end

    pyrAvgTrace = mean(pyrAll, 1, 'omitnan');
    intAvgTrace = mean(intAll, 1, 'omitnan');
end

function [xc, peakLagSec] = lagSweepTrialAverage(pyrAvgTrace, intAvgTrace, lags, pyrStart, pyrEnd, T, binSize)
    xc = nan(1, numel(lags));

    for iL = 1:numel(lags)
        L = lags(iL);

        intStart = pyrStart - L;
        intEnd   = pyrEnd   - L;

        if intStart < 1 || intEnd > T
            continue;
        end

        pseg = pyrAvgTrace(pyrStart:pyrEnd);
        iseg = intAvgTrace(intStart:intEnd);

        valid = ~isnan(pseg) & ~isnan(iseg);

        if nnz(valid) > 10
            xc(iL) = corr(pseg(valid)', iseg(valid)');
        end
    end

    [~, peakIdx] = max(xc);
    peakLagSec = lags(peakIdx) * binSize;
end

function xcZeroLag = zeroLagTrialAverage(pyrAvgTrace, intAvgTrace, pyrStart, pyrEnd)
    pseg = pyrAvgTrace(pyrStart:pyrEnd);
    iseg = intAvgTrace(pyrStart:pyrEnd);

    valid = ~isnan(pseg) & ~isnan(iseg);

    xcZeroLag = nan;
    if nnz(valid) > 10
        xcZeroLag = corr(pseg(valid)', iseg(valid)');
    end
end

function W = extractWindowsFromFR(frMat, neuralIdx1k, tAxis)
    nEvt  = numel(neuralIdx1k);
    nNeur = size(frMat,1);
    nTime = numel(tAxis);
    Ttot  = size(frMat,2);

    W = nan(nEvt, nNeur, nTime);

    for e = 1:nEvt
        t0 = neuralIdx1k(e);
        if isnan(t0)
            continue;
        end

        rng = t0 + tAxis;

        if any(rng < 1) || any(rng > Ttot)
            continue;
        end

        W(e,:,:) = frMat(:, rng);
    end
end

function [cortexIntFR, ok] = loadCortexIntFR(folderPath)
    ok = false;
    cortexIntFR = [];

    try
        frFile = fullfile(folderPath, 'NeuralFiringRates1msBins10msGauss.mat');
        if ~isfile(frFile), return; end
        load(frFile, 'cortexFRs', 'cortexInds');

        clsFile = fullfile('/home/asa7288/Transfer', 'AA_classifications.mat');
        if ~isfile(clsFile), return; end
        load(clsFile, 'classifications');

        if contains(folderPath,'D026'), matchRow = 1;
        elseif contains(folderPath,'D020'), matchRow = 2;
        elseif contains(folderPath,'D024'), matchRow = 3;
        elseif contains(folderPath,'D043'), matchRow = 4;
        else, return;
        end

        clsRow = classifications(matchRow, :);
        cortexIntFR = cortexFRs(clsRow{1,1}(cortexInds) == 1, :);
        ok = true;

    catch
        ok = false;
    end
end
