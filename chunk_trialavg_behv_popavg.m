% chunked TRIAL-AVERAGED population cross-correlation per behavior
% keeps:
%   - unpermuted lag-swept xc
%   - permuted lag-swept xc

% adds:
%   - trial averaging across behavior-specific events before cross-correlation
%   - shifted control (0-lag only, no lag sweep)

% shift behavior:
%   - first tries the precomputed shift for this jobInd
%   - if that does not yield a finite xcZeroLag for a behavior, samples a new shift on the fly
%   - keeps retrying until a finite xcZeroLag is obtained or maxRetries is reached
%   - saves the actual shift used

% IMPORTANT:
%   before rerunning shift jobs, delete old shift files in quest_runs
%   so the combine script does not mix old and new outputs

% jobInd meaning:
%   0 -> unpermuted real
%   1:100 -> permutation jobs
%   101:200 -> shifted-control jobs (requested shift #1..100)

clc;

%% ---------------- settings ----------------
baseDirs = {
    '/home/asa7288/Transfer/D026', ...
    '/home/asa7288/Transfer/D020', ...
    '/home/asa7288/Transfer/D024', ...
    '/home/asa7288/Transfer/D043', ...
    '/home/asa7288/Transfer/D050', ...
    '/home/asa7288/Transfer/D054'
};

labelsToUse = 'classifier';
labelTag = lower(labelsToUse);

rng(jobInd)

chunkHalf = 200;
maxLagSecs = 0.5;
doBaselineNorm = true;

switch lower(labelsToUse)
    case 'umap'
        behaviors = 1:7;
    case 'manual'
        behaviors = 1:10;
    case 'classifier'
        behaviors = 1:10;
    otherwise
        error('unknown labelsToUse: %s (use ''umap'', ''manual'', or ''classifier'')', labelsToUse);
end

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
    shiftInd = jobInd - 100;   % map 101:200 -> 1:100
else
    error('jobInd must be 0, 1:100, or 101:200');
end

for iDir = 1:numel(baseDirs)
    baseDir = baseDirs{iDir};
    fprintf('\nprocessing cortex session %d: %s (job %d, type=%s)\n', iDir, baseDir, jobInd, jobType);

    if jobType == "shift"
        S = load(fullfile(baseDir, 'EMG_Neural_AllChannels.mat'), 'pyrCxWinCell','intCxWinCell','tAxis', 'validTransitionsCell','validTransitionsNeurShiftedCell');

        syncData = load(fullfile(baseDir,'VideoSyncFrames.mat'), 'frameEMGSamples','frameNeuropixelSamples');
        frameEMGSamples = syncData.frameEMGSamples;
        frameNeuropixelSamples = syncData.frameNeuropixelSamples;

        emgData = load(fullfile(baseDir,'EMG1ms.mat'), 'downsampEMG');
        downsampEMG = emgData.downsampEMG;
        totalLength = size(downsampEMG, 2);

        minShift = 30000;  % 30 sec in ms
        maxRetries = 5000;
    else
        S = load(fullfile(baseDir, 'EMG_Neural_AllChannels.mat'), 'pyrCxWinCell','intCxWinCell','tAxis');
    end

    L = load(fullfile(baseDir, 'transitionCanonicalBehaviorLabels.mat'), 'regionLabelsPerTransition','manualLabelsPerTransition','classifierLabelsPerTransition');

    switch lower(labelsToUse)
        case 'umap'
            labelsPerTransition = L.regionLabelsPerTransition;
        case 'manual'
            labelsPerTransition = L.manualLabelsPerTransition;
        case 'classifier'
            labelsPerTransition = L.classifierLabelsPerTransition;
    end

    % EMG_Neural_AllChannels.mat now contains only good channels.
    % Use however many channels are actually saved for this session.
    channelsToUseThisSess = 1:numel(S.pyrCxWinCell);

    tAxis = S.tAxis(:)';
    binSize = 0.001;
    T = numel(tAxis);

    [~, cIdx] = min(abs(tAxis - 0));

    pyrStart = cIdx - chunkHalf;
    pyrEnd = cIdx + chunkHalf;

    reqMaxLagBins = round(maxLagSecs / binSize);
    maxLagLeft = pyrStart - 1;
    maxLagRight = T - pyrEnd;
    maxLagBins = min([reqMaxLagBins, maxLagLeft, maxLagRight]);
    lags = -maxLagBins:maxLagBins;

    for bIdx = 1:numel(behaviors)
        beh = behaviors(bIdx);
        fprintf('  behavior %d\n', beh);

        pyr_byCh = cell(1, numel(channelsToUseThisSess));
        int_byCh = cell(1, numel(channelsToUseThisSess));
        nTrialsThisBehavior = 0;

        if jobType == "shift"
            [cortexIntFR, okI] = loadCortexIntFR(baseDir);
            if ~okI
                warning('failed to load cortex interneuron FR for %s, skipping', baseDir);
                continue;
            end

            requestedShiftCol = shiftInd + 1;  % col1=real, col2:101=shift1:100
            usedShiftCol = NaN;
            usedShiftAmtMs = NaN;
            usedFallbackShift = false;
            foundValidShift = false;
            xcZeroLag = NaN;
            retryCount = 0;

            % ---------- first try the precomputed shift ----------
            [pyr_byCh_try, int_byCh_try, nTrials_try] = buildShiftedBehaviorSample( ...
                S, labelsPerTransition, beh, requestedShiftCol, ...
                cortexIntFR, tAxis, doBaselineNorm, false, [], [], [], [], [], channelsToUseThisSess);

            if nTrials_try > 0
                [pyrAvgTrace_try, intAvgTrace_try, nTrialsUsed_try] = ...
                    averageAcrossTrialsAndChannels(pyr_byCh_try, int_byCh_try);

                if nTrialsUsed_try > 0 && ~isempty(pyrAvgTrace_try) && ~isempty(intAvgTrace_try)
                    xcZeroLag_try = zeroLagTrialAverage( ...
                        pyrAvgTrace_try, intAvgTrace_try, pyrStart, pyrEnd);
                else
                    xcZeroLag_try = NaN;
                end

                if isfinite(xcZeroLag_try)
                    pyr_byCh = pyr_byCh_try;
                    int_byCh = int_byCh_try;
                    nTrialsThisBehavior = nTrials_try;
                    xcZeroLag = xcZeroLag_try;
                    usedShiftCol = requestedShiftCol;
                    usedFallbackShift = false;
                    foundValidShift = true;
                end
            end

            % ---------- fallback: sample new shifts on the fly until finite xcZeroLag ----------
            if ~foundValidShift
                for retryIdx = 1:maxRetries
                    retryCount = retryIdx;

                    randShift = randi([minShift, totalLength - minShift]);
                    if rand() > 0.5
                        randShift = -randShift;
                    end

                    [pyr_byCh_try, int_byCh_try, nTrials_try] = buildShiftedBehaviorSample( ...
                        S, labelsPerTransition, beh, NaN, ...
                        cortexIntFR, tAxis, doBaselineNorm, true, randShift, ...
                        frameEMGSamples, frameNeuropixelSamples, totalLength, baseDir, channelsToUseThisSess);

                    if nTrials_try < 1
                        continue;
                    end

                    [pyrAvgTrace_try, intAvgTrace_try, nTrialsUsed_try] = ...
                        averageAcrossTrialsAndChannels(pyr_byCh_try, int_byCh_try);

                    if nTrialsUsed_try > 0 && ~isempty(pyrAvgTrace_try) && ~isempty(intAvgTrace_try)
                        xcZeroLag_try = zeroLagTrialAverage( ...
                            pyrAvgTrace_try, intAvgTrace_try, pyrStart, pyrEnd);
                    else
                        xcZeroLag_try = NaN;
                    end

                    if isfinite(xcZeroLag_try)
                        pyr_byCh = pyr_byCh_try;
                        int_byCh = int_byCh_try;
                        nTrialsThisBehavior = nTrials_try;
                        xcZeroLag = xcZeroLag_try;
                        usedShiftCol = NaN;
                        usedShiftAmtMs = randShift;
                        usedFallbackShift = true;
                        foundValidShift = true;
                        break;
                    end
                end
            end

            if ~foundValidShift
                fprintf('    no valid shifted sample with finite xcZeroLag found for behavior %d after %d retries; skipping\n', ...
                    beh, maxRetries);
                continue;
            end

        else
            for ci = 1:numel(channelsToUseThisSess)
                ch = channelsToUseThisSess(ci);
                if ch > numel(labelsPerTransition) || isempty(labelsPerTransition{ch})
                    continue;
                end

                pyrWin = S.pyrCxWinCell{ch};
                intWin = S.intCxWinCell{ch};
                if isempty(pyrWin) || isempty(intWin)
                    continue;
                end

                labels = labelsPerTransition{ch}(:);
                if isempty(labels)
                    continue;
                end

                nEvtP = size(pyrWin,1);
                nEvtI = size(intWin,1);
                nLab  = numel(labels);
                nEvt  = min([nEvtP, nEvtI, nLab]);
                if nEvt < 1
                    continue;
                end

                labels = labels(1:nEvt);
                behMask = (labels == beh);
                if ~any(behMask)
                    continue;
                end

                evIdx = find(behMask);

                pyrWin_sel = pyrWin(evIdx,:,:);
                intWin_sel = intWin(evIdx,:,:);

                if jobType == "perm"
                    [nEvtP2, nPyr, ~] = size(pyrWin_sel);
                    [nEvtI2, nInt, ~] = size(intWin_sel);
                    nEvtUse = min(nEvtP2, nEvtI2);
                    if nPyr + nInt < 2 || nEvtUse < 1
                        continue;
                    end

                    pooledWin = cat(2, pyrWin_sel(1:nEvtUse,:,:), intWin_sel(1:nEvtUse,:,:));
                    totN = nPyr + nInt;
                    idx = randperm(totN);
                    idxP = idx(1:nPyr);
                    idxI = idx(nPyr+1:end);

                    pyrEvt = squeeze(mean(pooledWin(:, idxP, :), 2, 'omitnan'));
                    intEvt = squeeze(mean(pooledWin(:, idxI, :), 2, 'omitnan'));

                    if isvector(pyrEvt), pyrEvt = pyrEvt(:)'; end
                    if isvector(intEvt), intEvt = intEvt(:)'; end
                else
                    pyrEvt = meanEvt_keepEvents(pyrWin_sel);
                    intEvt = meanEvt_keepEvents(intWin_sel);
                end

                if doBaselineNorm
                    pyrEvt = subtractTrialBaseline(pyrEvt, tAxis, -500, -450);
                    intEvt = subtractTrialBaseline(intEvt, tAxis, -500, -450);
                end

                pyr_byCh{ci} = pyrEvt;
                int_byCh{ci} = intEvt;
                nTrialsThisBehavior = nTrialsThisBehavior + size(pyrEvt,1);
            end
        end

        if nTrialsThisBehavior == 0
            fprintf('    no trials for behavior %d; skipping\n', beh);
            continue;
        end

        tic

        % Trial-average behavior-specific event traces across all saved good channels.
        % Each row in pyr_byCh/int_byCh is already population-averaged across neurons
        % and baseline-subtracted on a per-event basis.
        [pyrAvgTrace, intAvgTrace, nTrialsUsed] = ...
            averageAcrossTrialsAndChannels(pyr_byCh, int_byCh);

        if isempty(pyrAvgTrace) || isempty(intAvgTrace) || nTrialsUsed < 1
            fprintf('    no usable trial-averaged traces for behavior %d; skipping\n', beh);
            continue;
        end

        if jobType == "shift"
            xcZeroLag = zeroLagTrialAverage( ...
                pyrAvgTrace, intAvgTrace, pyrStart, pyrEnd);

            runtimeSec = toc;
            fprintf('    saving shift file: beh=%d, shiftInd=%d, fallback=%d, xcZeroLag=%.4f, nTrials=%d, retries=%d\n', ...
                beh, shiftInd, usedFallbackShift, xcZeroLag, nTrialsUsed, retryCount);
            disp(['shifted trial-averaged control done, time: ' num2str(runtimeSec) ' s'])
        else
            [xc, peakLagSec] = lagSweepTrialAverage( ...
                pyrAvgTrace, intAvgTrace, lags, pyrStart, pyrEnd, T, binSize);

            runtimeSec = toc;
            disp(['trial-averaged cross-correlation analysis done, time: ' num2str(runtimeSec) ' s'])
        end

        outDir = fullfile(baseDir, 'quest_runs');
        if ~exist(outDir, 'dir'), mkdir(outDir); end

        if jobType == "real"
            outFile = fullfile(outDir, ...
                sprintf('concatCrossCorrPerCanonicalBehavior_%s_unperm_beh%02d.mat', labelTag, beh));
            save(outFile, ...
                'lags','binSize','xc','peakLagSec','chunkHalf', ...
                'doBaselineNorm','jobInd','permInd','shiftInd','jobType', ...
                'baseDir','runtimeSec','beh','nTrialsThisBehavior','nTrialsUsed','labelsToUse', ...
                'pyrAvgTrace','intAvgTrace','channelsToUseThisSess');

        elseif jobType == "perm"
            outFile = fullfile(outDir, ...
                sprintf('concatCrossCorrPerCanonicalBehavior_%s_perm_%03d_beh%02d.mat', labelTag, permInd, beh));
            save(outFile, ...
                'lags','binSize','xc','peakLagSec','chunkHalf', ...
                'doBaselineNorm','jobInd','permInd','shiftInd','jobType', ...
                'baseDir','runtimeSec','beh','nTrialsThisBehavior','nTrialsUsed','labelsToUse', ...
                'pyrAvgTrace','intAvgTrace','channelsToUseThisSess');

        else
            outFile = fullfile(outDir, ...
                sprintf('concatCrossCorrPerCanonicalBehavior_%s_shift_%03d_zerolag_beh%02d.mat', labelTag, shiftInd, beh));
            save(outFile, ...
                'binSize','xcZeroLag','chunkHalf', ...
                'doBaselineNorm','jobInd','permInd','shiftInd','jobType', ...
                'baseDir','runtimeSec','beh','nTrialsThisBehavior','nTrialsUsed','labelsToUse', ...
                'requestedShiftCol','usedShiftCol','usedShiftAmtMs','usedFallbackShift','maxRetries','retryCount', ...
                'pyrAvgTrace','intAvgTrace','channelsToUseThisSess');
        end

        fprintf('saved: %s\n', outFile);
    end
end

%% ===================== helpers =====================

function [pyr_byCh, int_byCh, nTrialsThisBehavior] = buildShiftedBehaviorSample( ...
    S, labelsPerTransition, beh, shiftCol, cortexIntFR, tAxis, doBaselineNorm, ...
    useFallbackShift, randShift, frameEMGSamples, frameNeuropixelSamples, totalLength, baseDir, channelsToUseThisSess)

    pyr_byCh = cell(1, numel(channelsToUseThisSess));
    int_byCh = cell(1, numel(channelsToUseThisSess));
    nTrialsThisBehavior = 0;

    for ci = 1:numel(channelsToUseThisSess)
        ch = channelsToUseThisSess(ci);
        if ch > numel(labelsPerTransition) || isempty(labelsPerTransition{ch})
            continue;
        end

        pyrWin = S.pyrCxWinCell{ch};
        intWin = S.intCxWinCell{ch};

        if isempty(pyrWin) || isempty(intWin)
            continue;
        end

        labels = labelsPerTransition{ch}(:);
        if isempty(labels)
            continue;
        end

        nEvtP = size(pyrWin,1);
        nEvtI = size(intWin,1);
        nLab  = numel(labels);
        nEvt  = min([nEvtP, nEvtI, nLab]);

        if nEvt < 1
            continue;
        end

        labels = labels(1:nEvt);
        behMask = (labels == beh);
        if ~any(behMask)
            continue;
        end

        evIdx = find(behMask);
        pyrWin_sel = pyrWin(evIdx,:,:);

        if ~useFallbackShift
            idx1k_all = S.validTransitionsNeurShiftedCell{ch, shiftCol};
            if isempty(idx1k_all)
                continue;
            end

            idx1k_all = idx1k_all(1:nEvt);
            idx1k = idx1k_all(evIdx);

        else
            validTransitions = S.validTransitionsCell{ch};
            if isempty(validTransitions)
                continue;
            end

            validTransitions = validTransitions(1:nEvt);
            transitionsThisBeh = validTransitions(evIdx);

            shifted = transitionsThisBeh + randShift;
            shifted(shifted > totalLength) = shifted(shifted > totalLength) - (totalLength + 1);
            shifted(shifted < 1) = shifted(shifted < 1) + totalLength;

            currentDir = pwd;
            cd(baseDir);
            neuralIndices = NeurEMGSync(shifted * 20, frameEMGSamples, frameNeuropixelSamples, 'EMG');
            cd(currentDir);

            idx1k = round(neuralIndices / 30);
        end

        validEvtMask = ~isnan(idx1k);
        idx1k = idx1k(validEvtMask);
        pyrWin_sel = pyrWin_sel(validEvtMask,:,:);

        if isempty(idx1k) || isempty(pyrWin_sel)
            continue;
        end

        pyrEvt = meanEvt_keepEvents(pyrWin_sel);
        intWinShift = extractWindowsFromFR(cortexIntFR, idx1k, tAxis);
        intEvt = meanEvt_keepEvents(intWinShift);

        nEvtUse = min(size(pyrEvt,1), size(intEvt,1));
        if nEvtUse < 1
            continue;
        end

        pyrEvt = pyrEvt(1:nEvtUse,:);
        intEvt = intEvt(1:nEvtUse,:);

        if doBaselineNorm
            pyrEvt = subtractTrialBaseline(pyrEvt, tAxis, -500, -450);
            intEvt = subtractTrialBaseline(intEvt, tAxis, -500, -450);
        end

        pyr_byCh{ci} = pyrEvt;
        int_byCh{ci} = intEvt;
        nTrialsThisBehavior = nTrialsThisBehavior + size(pyrEvt,1);
    end
end

function Eout = subtractTrialBaseline(Ein, tAxisSec, tStartSec, tEndSec)
    [~, iStart] = min(abs(tAxisSec - tStartSec));
    [~, iEnd]   = min(abs(tAxisSec - tEndSec));
    if iEnd < iStart
        tmp = iStart; iStart = iEnd; iEnd = tmp;
    end
    baselines = mean(Ein(:, iStart:iEnd), 2, 'omitnan');
    Eout = Ein - baselines;
end

function E = meanEvt_keepEvents(M)
    if isempty(M)
        E = [];
        return;
    end
    E = mean(M, 2, 'omitnan');
    E = permute(E, [1 3 2]);
    E = reshape(E, size(E,1), size(E,2));
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
        intEnd   = pyrEnd - L;

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
        if isnan(t0), continue; end
        rng = t0 + tAxis;
        if any(rng < 1) || any(rng > Ttot), continue; end
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
        elseif contains(folderPath,'D050'), matchRow = 5;
        elseif contains(folderPath,'D054'), matchRow = 6;
        else, return; end

        clsRow = classifications(matchRow, :);
        cortexIntFR = cortexFRs(clsRow{1,1}(cortexInds) == 1, :);
        ok = true;

    catch
        ok = false;
    end
end
