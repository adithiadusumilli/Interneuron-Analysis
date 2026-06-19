% pairwise_chunked_trialavg_allpairs_xcorr_quest.m

% Trial-averaged pairwise xcorr with proper shift null

% jobInd = 0 -> real all-pairs trial-averaged xcorr

% jobInd = 1:100 -> one null matrix per shift job, exactly matching the old single-trial shift pipeline but with trial averaging:

% OLD: corr(allUnshift(trials,i,:), allShift(trials,j,:))
% concatenated across trials, using physically shifted event timestamps

% NEW: corr(mean(allUnshift(:,i,:)), mean(allShift(:,j,:)))
% same physically shifted event timestamps from validTransitionsNeurShiftedCell, same FR extraction via extractWindowsFromFR, j avg across trials before corr

% shift is identical to the old pipeline, only corr step changes from single-trial concatenation to trial-avg traces

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

assert(exist('jobInd','var') == 1, 'define jobInd before running');
assert(jobInd >= 0 && jobInd <= 100, 'jobInd must be 0:100');

isShiftedJob = (jobInd > 0);

rng(jobInd);

for iDir = 1:numel(baseDirs)

    baseDir = baseDirs{iDir};
    fprintf('\nprocessing %s | sessInd=%d | jobInd=%d\n', baseDir, sessInd, jobInd);

    %% load data — same fields as old pipeline
    S = load(fullfile(baseDir, 'EMG_Neural_AllChannels.mat'), ...
        'pyrCxWinCell', 'intCxWinCell', 'tAxis', ...
        'validTransitionsNeurCell', 'validTransitionsNeurShiftedCell');

    if ~isfield(S,'tAxis') || ~isfield(S,'pyrCxWinCell') || ~isfield(S,'intCxWinCell')
        fprintf('missing required vars in %s, skipping\n', baseDir);
        continue;
    end

    if isShiftedJob
        if ~isfield(S, 'validTransitionsNeurShiftedCell')
            fprintf('missing validTransitionsNeurShiftedCell in %s, skipping\n', baseDir);
            continue;
        end
        if size(S.validTransitionsNeurShiftedCell, 2) < 101
            error('expected validTransitionsNeurShiftedCell to have 101 columns (real + 100 shifts)');
        end
    end

    tAxis = S.tAxis(:)';
    binSize = 0.001;
    T = numel(tAxis);

    [~, cIdx] = min(abs(tAxis - 0));
    chunkStart = cIdx - chunkHalf;
    chunkEnd = cIdx + chunkHalf;

    if chunkStart < 1 || chunkEnd > T
        error('chunkHalf=%d does not fit in window length=%d', chunkHalf, T);
    end

    reqMaxLagBins = round(maxLagSecs / binSize);
    maxLagLeft = chunkStart - 1;
    maxLagRight = T - chunkEnd;
    maxLagBins = min([reqMaxLagBins, maxLagLeft, maxLagRight]);
    lags = -maxLagBins:maxLagBins;

    %% build unshifted pooled windows across channels (same as old pipeline)
    pyrAll_unshift = [];
    intAll_unshift = [];
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
        pyrWin = pyrWin(1:nEvt, :, :);
        intWin = intWin(1:nEvt, :, :);

        if isnan(nPyr_ref), nPyr_ref = size(pyrWin,2); end
        if isnan(nInt_ref), nInt_ref = size(intWin,2); end

        nPyr_ref = min(nPyr_ref, size(pyrWin,2));
        nInt_ref = min(nInt_ref, size(intWin,2));

        pyrWin = pyrWin(:, 1:nPyr_ref, :);
        intWin = intWin(:, 1:nInt_ref, :);

        if doBaselineNorm
            pyrWin = subtractTrialBaseline3d(pyrWin, tAxis, -500, -450);
            intWin = subtractTrialBaseline3d(intWin, tAxis, -500, -450);
        end

        pyrAll_unshift = cat(1, pyrAll_unshift, pyrWin);
        intAll_unshift = cat(1, intAll_unshift, intWin);
    end

    if isempty(pyrAll_unshift) || isempty(intAll_unshift)
        fprintf('no valid channels found, skipping\n');
        continue;
    end

    nEvt_un        = min(size(pyrAll_unshift,1), size(intAll_unshift,1));
    pyrAll_unshift = pyrAll_unshift(1:nEvt_un, :, :);
    intAll_unshift = intAll_unshift(1:nEvt_un, :, :);

    % ordering [pyr, int] to match old pipeline
    allAll_unshift = cat(2, pyrAll_unshift, intAll_unshift);  % [nEvt x nAll x nTime]
    nAll           = size(allAll_unshift, 2);

    neuronType = [zeros(1,nPyr_ref), ones(1,nInt_ref)];
    pyrIdx     = 1:nPyr_ref;
    intIdx     = (nPyr_ref+1):nAll;

    outDir = fullfile(baseDir, 'quest_runs');
    if ~exist(outDir,'dir'), mkdir(outDir); end

    tic;

    if ~isShiftedJob
        %% ---- jobInd=0: trial-average then pairwise xcorr ----

        % trial-average the unshifted windows
        % ordering [pyr, int] preserved
        pyrAvg = squeeze(mean(pyrAll_unshift, 1, 'omitnan'));  % [nPyr x nTime]
        intAvg = squeeze(mean(intAll_unshift, 1, 'omitnan'));  % [nInt x nTime]
        allAvg = cat(1, pyrAvg, intAvg);                       % [nAll x nTime]

        fprintf('real job: nEvt=%d | nAll=%d (nPyr=%d, nInt=%d)\n', ...
            nEvt_un, nAll, nPyr_ref, nInt_ref);

        [xcMat_all, peakCorrMat_all, peakLagSecMat_all] = ...
            lagSweepPairwise_trialavg_allpairs( ...
                allAvg, lags, chunkStart, chunkEnd, T, binSize);

        runtimeSec = toc;
        fprintf('real all-pairs done: %.1f s\n', runtimeSec);

        outFile = fullfile(outDir, ...
            sprintf('pairwiseChunkedTrialAvgXCorr_ALLPAIRS_unshifted_sess%02d_fullxc.mat', sessInd));

        save(outFile, ...
            'lags', 'binSize', 'tAxis', ...
            'xcMat_all', 'peakCorrMat_all', 'peakLagSecMat_all', ...
            'chunkHalf', 'channelsToUse', 'doBaselineNorm', ...
            'jobInd', 'baseDir', 'sessInd', ...
            'nPyr_ref', 'nInt_ref', 'nAll', 'neuronType', 'pyrIdx', 'intIdx', ...
            'runtimeSec', '-v7.3');

    else
        %% ---- jobInd=1:100: extract shifted windows (same as old pipeline)
        %%                     then trial-average, then xcorr ----
        %
        % Identical to old pipeline up to the averaging step:
        %   1. load cortex FR matrices (same files, same paths)
        %   2. use validTransitionsNeurShiftedCell{ch, shiftCol} (same shifted timestamps)
        %   3. extract windows via extractWindowsFromFR (same function)
        %   4. baseline normalize (same)
        %   NEW: trial-average the shifted windows before correlating
        %   NEW: xcorr on averaged traces instead of concatenated trials

        [cortexIntFR, okI] = loadCortexIntFR(baseDir);
        [cortexPyrFR, okP] = loadCortexPyrFR(baseDir);

        if ~okI || ~okP
            fprintf('failed to load cortex FR for %s, skipping\n', baseDir);
            continue;
        end

        shiftCol = jobInd + 1;

        pyrAll_shift = [];
        intAll_shift = [];

        for ci = 1:numel(channelsToUse)
            ch = channelsToUse(ci);

            idx1k = S.validTransitionsNeurShiftedCell{ch, shiftCol};
            if isempty(idx1k)
                continue;
            end

            pyrWinShift = extractWindowsFromFR(cortexPyrFR, idx1k, tAxis);
            intWinShift = extractWindowsFromFR(cortexIntFR, idx1k, tAxis);

            % match neuron counts to unshifted refs (same as old pipeline)
            if size(pyrWinShift,2) < nPyr_ref
                nPyr_ref       = size(pyrWinShift,2);
                allAll_unshift = allAll_unshift(:, [1:nPyr_ref, (nAll-nInt_ref+1):nAll], :);
                pyrIdx         = 1:nPyr_ref;
                nAll           = nPyr_ref + nInt_ref;
            end
            if size(intWinShift,2) < nInt_ref
                nInt_ref       = size(intWinShift,2);
                allAll_unshift = allAll_unshift(:, 1:(nPyr_ref+nInt_ref), :);
                intIdx         = (nPyr_ref+1):(nPyr_ref+nInt_ref);
                nAll           = nPyr_ref + nInt_ref;
            end

            pyrWinShift = pyrWinShift(:, 1:nPyr_ref, :);
            intWinShift = intWinShift(:, 1:nInt_ref, :);

            if doBaselineNorm
                pyrWinShift = subtractTrialBaseline3d(pyrWinShift, tAxis, -500, -450);
                intWinShift = subtractTrialBaseline3d(intWinShift, tAxis, -500, -450);
            end

            pyrAll_shift = cat(1, pyrAll_shift, pyrWinShift);
            intAll_shift = cat(1, intAll_shift, intWinShift);
        end

        if isempty(pyrAll_shift) || isempty(intAll_shift)
            fprintf('no shifted events found for jobInd=%d, skipping\n', jobInd);
            continue;
        end

        % match event counts between unshifted and shifted (same as old pipeline)
        nEvt_sh      = min(size(pyrAll_shift,1), size(intAll_shift,1));
        pyrAll_shift = pyrAll_shift(1:nEvt_sh, :, :);
        intAll_shift = intAll_shift(1:nEvt_sh, :, :);

        allAll_shift = cat(2, pyrAll_shift, intAll_shift);   % [nEvt x nAll x nTime]

        nEvt           = min(size(allAll_unshift,1), size(allAll_shift,1));
        allAll_unshift = allAll_unshift(1:nEvt, :, :);
        allAll_shift   = allAll_shift(1:nEvt, :, :);
        nAll           = size(allAll_unshift, 2);

        fprintf('shift job: nEvt=%d | nAll=%d (nPyr=%d, nInt=%d)\n', ...
            nEvt, nAll, nPyr_ref, nInt_ref);

        %% trial-average both unshifted and shifted windows
        % NEW vs old pipeline: average here before correlating
        allAvgUnshift = squeeze(mean(allAll_unshift, 1, 'omitnan'));  % [nAll x nTime]
        allAvgShift   = squeeze(mean(allAll_shift,   1, 'omitnan'));  % [nAll x nTime]

        %% compute null: corr(unshifted_avg_i, shifted_avg_j) for each pair (i<j)
        % directly mirrors old: corr(allUnshift(trials,i,:), allShift(trials,j,:))
        % just averaged first
        nullCorrMat_all = pairwise_trialavg_unshift_vs_shift( ...
            allAvgUnshift, allAvgShift, lags, chunkStart, chunkEnd, T, binSize);

        runtimeSec = toc;
        fprintf('shift null done: %.1f s\n', runtimeSec);

        outFile = fullfile(outDir, ...
            sprintf('pairwiseChunkedTrialAvgXCorr_ALLPAIRS_sess%02d_shiftNull_%03d.mat', sessInd, jobInd));

        save(outFile, ...
            'nullCorrMat_all', ...
            'lags', 'binSize', 'tAxis', ...
            'chunkHalf', 'channelsToUse', 'doBaselineNorm', ...
            'jobInd', 'baseDir', 'sessInd', ...
            'nPyr_ref', 'nInt_ref', 'nAll', 'neuronType', 'pyrIdx', 'intIdx', ...
            'runtimeSec', '-v7.3');
    end

    fprintf('saved: %s\n', outFile);
end

%% ===================== helpers =====================

function Wout = subtractTrialBaseline3d(Win, tAxis, tStart, tEnd)
    [~, iStart] = min(abs(tAxis - tStart));
    [~, iEnd]   = min(abs(tAxis - tEnd));
    if iEnd < iStart, tmp = iStart; iStart = iEnd; iEnd = tmp; end
    base = mean(Win(:,:,iStart:iEnd), 3, 'omitnan');
    Wout = Win - base;
end

function [xcMat, peakCorrMat, peakLagSecMat] = lagSweepPairwise_trialavg_allpairs( ...
    allAvg, lags, chunkStart, chunkEnd, T, binSize)
% Pairwise lag sweep on trial-averaged traces [nAll x nTime].
% Upper triangle only (j>i).

    nAll = size(allAvg, 1);
    nL   = numel(lags);

    xcMat         = nan(nAll, nAll, nL);
    peakCorrMat   = nan(nAll, nAll);
    peakLagSecMat = nan(nAll, nAll);

    for i = 1:nAll
        for j = (i+1):nAll
            [xc, peakLagSec, peakCorr] = onePairLagSweep( ...
                allAvg(i,:), allAvg(j,:), lags, chunkStart, chunkEnd, T, binSize);

            xcMat(i,j,:)       = xc;
            peakCorrMat(i,j)   = peakCorr;
            peakLagSecMat(i,j) = peakLagSec;
        end
    end
end

function nullCorrMat = pairwise_trialavg_unshift_vs_shift( ...
    allAvgUnshift, allAvgShift, lags, chunkStart, chunkEnd, T, binSize)
% Null correlation matrix using trial-averaged traces.
% For each pair (i<j): corr(unshifted_avg_i, shifted_avg_j) at peak lag.
% Mirrors old: corr(allUnshift(trials,i,:), allShift(trials,j,:))
% Mirror to lower triangle to match old pipeline symmetric structure.

    nAll        = size(allAvgUnshift, 1);
    nullCorrMat = nan(nAll, nAll);

    for i = 1:nAll
        for j = (i+1):nAll
            [~, ~, peakCorr] = onePairLagSweep( ...
                allAvgUnshift(i,:), allAvgShift(j,:), ...
                lags, chunkStart, chunkEnd, T, binSize);

            nullCorrMat(i,j) = peakCorr;
            nullCorrMat(j,i) = peakCorr;
        end
    end
end

function [xc, peakLagSec, peakCorr] = onePairLagSweep( ...
    aTrace, bTrace, lags, chunkStart, chunkEnd, T, binSize)

    xc = nan(1, numel(lags));

    for iL = 1:numel(lags)
        L      = lags(iL);
        bStart = chunkStart - L;
        bEnd   = chunkEnd   - L;

        if bStart < 1 || bEnd > T, continue; end

        aseg  = aTrace(chunkStart:chunkEnd);
        bseg  = bTrace(bStart:bEnd);
        valid = ~isnan(aseg) & ~isnan(bseg);

        if nnz(valid) > 10
            xc(iL) = corr(aseg(valid)', bseg(valid)');
        end
    end

    [peakCorr, peakIdx] = max(xc);
    peakLagSec          = lags(peakIdx) * binSize;
end

function W = extractWindowsFromFR(frMat, neuralIdx1k, tAxis)
    nEvt  = numel(neuralIdx1k);
    nNeur = size(frMat, 1);
    nTime = numel(tAxis);
    Ttot  = size(frMat, 2);

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
    ok = false; cortexIntFR = [];
    try
        frFile = fullfile(folderPath, 'NeuralFiringRates1msBins10msGauss.mat');
        load(frFile, 'cortexFRs', 'cortexInds');

        clsFile = fullfile('/home/asa7288/Transfer', 'AA_classifications.mat');
        load(clsFile, 'classifications');

        if contains(folderPath,'D026'),      matchRow = 1;
        elseif contains(folderPath,'D020'), matchRow = 2;
        elseif contains(folderPath,'D024'), matchRow = 3;
        elseif contains(folderPath,'D043'), matchRow = 4;
        else, error('unknown session folder'); end

        cls          = classifications{matchRow,1}(cortexInds);
        cortexIntFR  = cortexFRs(cls == 1, :);
        ok           = true;
    catch ME
        warning(ME.message);
    end
end

function [cortexPyrFR, ok] = loadCortexPyrFR(folderPath)
    ok = false; cortexPyrFR = [];
    try
        frFile = fullfile(folderPath, 'NeuralFiringRates1msBins10msGauss.mat');
        load(frFile, 'cortexFRs', 'cortexInds');

        clsFile = fullfile('/home/asa7288/Transfer', 'AA_classifications.mat');
        load(clsFile, 'classifications');

        if contains(folderPath,'D026'),      matchRow = 1;
        elseif contains(folderPath,'D020'), matchRow = 2;
        elseif contains(folderPath,'D024'), matchRow = 3;
        elseif contains(folderPath,'D043'), matchRow = 4;
        else, error('unknown session folder'); end

        cls          = classifications{matchRow,1}(cortexInds);
        cortexPyrFR  = cortexFRs(cls == 0, :);
        ok           = true;
    catch ME
        warning(ME.message);
    end
end
