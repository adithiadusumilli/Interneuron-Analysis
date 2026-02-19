% pairwise_chunked_allpairs_xcorr_with_shifting_quest.m

% chunked ALL-PAIRS pairwise xcorr for cortex:
%   - int-int, pyr-pyr, int-pyr
% no duplicates:
%   - REAL (jobInd=0): computes only upper triangle (i<j) and saves full xc traces + peaks
%   - SHIFT (jobInd=1..100): computes ONE null matrix per shift where each pair is:
%         (unshifted i) vs (shifted j) for i<j
%     then mirrors to keep symmetric

clearvars; clc;

%% ---------------- settings ----------------
baseDirs = {
    '/home/asa7288/Transfer/D026', ...
    '/home/asa7288/Transfer/D020', ...
    '/home/asa7288/Transfer/D024', ...
    '/home/asa7288/Transfer/D043'
};

sessInd = 1; 
baseDirs = baseDirs(sessInd); % reduce to 1 animal

channelsToUse  = 1:4;
chunkHalf = 200;   % 401 samples
maxLagSecs = 0.5;  % only used for unshifted
doBaselineNorm = true;

%% ---------------- main loop ----------------
assert(exist('jobInd','var')==1, 'define jobInd before running (0=real, 1..100=shift)');
assert(jobInd >= 0 && jobInd <= 100, 'jobInd must be 0..100');
isShiftedJob = (jobInd > 0);

for iDir = 1:numel(baseDirs)
    baseDir = baseDirs{iDir};
    fprintf('\nprocessing cortex session %d: %s (job %d)\n', iDir, baseDir, jobInd);

    % -------- load windowed data + indices --------
    S = load(fullfile(baseDir, 'EMG_Neural_AllChannels.mat'), 'pyrCxWinCell', 'intCxWinCell', 'tAxis', 'validTransitionsNeurCell', 'validTransitionsNeurShiftedCell');

    if ~isfield(S,'tAxis') || ~isfield(S,'pyrCxWinCell') || ~isfield(S,'intCxWinCell')
        fprintf('missing required vars in %s, skipping\n', baseDir);
        continue;
    end

    if isShiftedJob
        if ~isfield(S,'validTransitionsNeurShiftedCell')
            fprintf('missing validTransitionsNeurShiftedCell in %s, skipping shifted job\n', baseDir);
            continue;
        end
        if size(S.validTransitionsNeurShiftedCell,2) < 101
            error('expected validTransitionsNeurShiftedCell to have 101 columns (real + 100 shifts).');
        end
    end

    % -------- window geometry --------
    tAxis = S.tAxis(:)'; % ms offsets, length ~1001
    binSize = 0.001; % seconds
    T = numel(tAxis);

    [~, cIdx] = min(abs(tAxis - 0));
    chunkStart = cIdx - chunkHalf;
    chunkEnd = cIdx + chunkHalf;

    if chunkStart < 1 || chunkEnd > T
        error('chunkHalf=%d does not fit in window length=%d.', chunkHalf, T);
    end

    % lag sweep range (unshifted only)
    reqMaxLagBins = round(maxLagSecs / binSize);
    maxLagLeft = chunkStart - 1;
    maxLagRight = T - chunkEnd;
    maxLagBins = min([reqMaxLagBins, maxLagLeft, maxLagRight]);

    lags = -maxLagBins:maxLagBins;  % bins
    nL = numel(lags); %#ok<NASGU>

    % ============================================================
    % build UNHIFTED pooled windows for pyr + int across channels
    % ============================================================
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

        % match events across types per channel
        nEvt = min(size(pyrWin,1), size(intWin,1));
        pyrWin = pyrWin(1:nEvt,:,:);
        intWin = intWin(1:nEvt,:,:);

        % enforce consistent neuron counts across channels (separately for pyr/int)
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

        pyrAll_unshift = cat(1, pyrAll_unshift, pyrWin);
        intAll_unshift = cat(1, intAll_unshift, intWin);
    end

    if isempty(pyrAll_unshift) || isempty(intAll_unshift)
        fprintf('no valid channels found, skipping\n');
        continue;
    end

    % match event counts in unshifted
    nEvt_un = min(size(pyrAll_unshift,1), size(intAll_unshift,1));
    pyrAll_unshift = pyrAll_unshift(1:nEvt_un,:,:);
    intAll_unshift = intAll_unshift(1:nEvt_un,:,:);

    % combine unshifted into ALL ordering: [pyr, int]
    allAll_unshift = cat(2, pyrAll_unshift, intAll_unshift);  % events x nAll x time
    nAll = size(allAll_unshift,2);

    neuronType = [zeros(1,nPyr_ref), ones(1,nInt_ref)];  % 0=pyr, 1=int
    pyrIdx = 1:nPyr_ref;
    intIdx = (nPyr_ref+1):nAll;

    % ============================================================
    % SHIFT jobs: build SHIFTED pooled windows from FR, same ordering
    % ============================================================
    if isShiftedJob
        [cortexIntFR, okI] = loadCortexIntFR(baseDir);
        [cortexPyrFR, okP] = loadCortexPyrFR(baseDir);
        if ~okI || ~okP
            fprintf('failed to load cortex pyr/int FR for %s, skipping\n', baseDir);
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

            % match neuron counts to unshifted refs
            if size(pyrWinShift,2) < nPyr_ref
                nPyr_ref = size(pyrWinShift,2);
                allAll_unshift = allAll_unshift(:, [1:nPyr_ref, (size(allAll_unshift,2)-nInt_ref+1):size(allAll_unshift,2)], :);
                pyrIdx = 1:nPyr_ref;
            end
            if size(intWinShift,2) < nInt_ref
                nInt_ref = size(intWinShift,2);
                allAll_unshift = allAll_unshift(:, 1:(nPyr_ref+nInt_ref), :);
                intIdx = (nPyr_ref+1):(nPyr_ref+nInt_ref);
            end

            pyrWinShift = pyrWinShift(:,1:nPyr_ref,:);
            intWinShift = intWinShift(:,1:nInt_ref,:);

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

        % match event counts in shifted
        nEvt_sh = min(size(pyrAll_shift,1), size(intAll_shift,1));
        pyrAll_shift = pyrAll_shift(1:nEvt_sh,:,:);
        intAll_shift = intAll_shift(1:nEvt_sh,:,:);

        allAll_shift = cat(2, pyrAll_shift, intAll_shift); % events x nAll x time

        % final safety: match event counts between unshifted and shifted
        nEvt = min(size(allAll_unshift,1), size(allAll_shift,1));
        allAll_unshift = allAll_unshift(1:nEvt,:,:);
        allAll_shift   = allAll_shift(1:nEvt,:,:);

        nAll = size(allAll_unshift,2);
        fprintf('shift job: nEvt=%d | nAll=%d (nPyr=%d, nInt=%d)\n', nEvt, nAll, nPyr_ref, nInt_ref);

    else
        nEvt = size(allAll_unshift,1);
        fprintf('real job: nEvt=%d | nAll=%d (nPyr=%d, nInt=%d)\n', nEvt, nAll, nPyr_ref, nInt_ref);
    end

    % ============================================================
    % compute metrics
    % ============================================================
    tic
    if ~isShiftedJob
        [xcMat_all, peakCorrMat_all, peakLagSecMat_all] = ...
            lagSweepPairwise_fullxc_chunked_allpairs(allAll_unshift, lags, chunkStart, chunkEnd, T, binSize);

        runtimeSec = toc;
        disp(['real all-pairs done, time: ' num2str(runtimeSec) ' s'])

    else
        % KEY CHANGE:
        % null for each pair (i<j): (unshifted i) vs (shifted j)
        nullCorrMat_all = pairwise_zeroLag_chunked_allpairs_unshift_vs_shift(allAll_unshift, allAll_shift, chunkStart, chunkEnd);

        runtimeSec = toc;
        disp(['shift all-pairs done, time: ' num2str(runtimeSec) ' s'])
    end

    % ============================================================
    % save
    % ============================================================
    outDir = fullfile(baseDir, 'quest_runs');
    if ~exist(outDir,'dir'), mkdir(outDir); end

    if ~isShiftedJob
        outFile = fullfile(outDir, sprintf('pairwiseChunkedXCorr_ALLPAIRS_unshifted_sess%02d_fullxc.mat', sessInd));
        save(outFile, 'lags','binSize','tAxis', 'xcMat_all','peakCorrMat_all','peakLagSecMat_all', 'chunkHalf','channelsToUse','doBaselineNorm', 'jobInd','baseDir','sessInd', 'nPyr_ref','nInt_ref','nAll','neuronType','pyrIdx','intIdx', 'runtimeSec','-v7.3');
    else
        outFile = fullfile(outDir, sprintf('pairwiseChunkedXCorr_ALLPAIRS_sess%02d_shift_%03d_zerolag.mat', sessInd, jobInd));
        save(outFile, 'binSize','tAxis', 'nullCorrMat_all', 'chunkHalf','channelsToUse','doBaselineNorm', 'jobInd','baseDir','sessInd', 'nPyr_ref','nInt_ref','nAll','neuronType','pyrIdx','intIdx', 'runtimeSec', '-v7.3');
    end

    fprintf('saved: %s\n', outFile);
end

%% ===================== helpers =====================

function Wout = subtractTrialBaseline3d(Win, tAxis, tStart, tEnd)
    [~, iStart] = min(abs(tAxis - tStart));
    [~, iEnd]   = min(abs(tAxis - tEnd));
    if iEnd < iStart, tmp=iStart; iStart=iEnd; iEnd=tmp; end
    base = mean(Win(:,:,iStart:iEnd), 3, 'omitnan');
    Wout = Win - base;
end

function [xcMat, peakCorrMat, peakLagSecMat] = lagSweepPairwise_fullxc_chunked_allpairs(allAll, lags, chunkStart, chunkEnd, T, binSize)
% allAll: events x nAll x time
% computes only upper triangle (j>i), lower triangle stays NaN
    nEvt = size(allAll,1);
    nAll = size(allAll,2);
    nL = numel(lags);
    segLen = chunkEnd - chunkStart + 1;

    xcMat = nan(nAll, nAll, nL);
    peakCorrMat = nan(nAll, nAll);
    peakLagSecMat = nan(nAll, nAll);

    for i = 1:nAll
        for j = (i+1):nAll

            bestR = -inf;
            bestLagSec = nan;

            for iL = 1:nL
                L = lags(iL);

                jStart = chunkStart - L;
                jEnd   = chunkEnd   - L;

                if jStart < 1 || jEnd > T
                    continue;
                end

                Avec = nan(1, nEvt * segLen);
                Bvec = nan(1, nEvt * segLen);

                w = 0;
                for e = 1:nEvt
                    aseg = squeeze(allAll(e, i, chunkStart:chunkEnd));
                    bseg = squeeze(allAll(e, j, jStart:jEnd));

                    idx = (w+1):(w+segLen);
                    Avec(idx) = aseg(:)';
                    Bvec(idx) = bseg(:)';
                    w = w + segLen;
                end

                valid = ~isnan(Avec) & ~isnan(Bvec);
                if nnz(valid) > 10
                    r = corr(Avec(valid)', Bvec(valid)');
                    xcMat(i,j,iL) = r;

                    if ~isnan(r) && r > bestR
                        bestR = r;
                        bestLagSec = L * binSize;
                    end
                end
            end

            if isfinite(bestR)
                peakCorrMat(i,j) = bestR;
                peakLagSecMat(i,j) = bestLagSec;
            end
        end
    end
end

function nullCorrMat = pairwise_zeroLag_chunked_allpairs_unshift_vs_shift(allUnshift, allShift, chunkStart, chunkEnd)
% KEY NULL:
% for each i<j:
%   use unshifted neuron i chunk vs shifted neuron j chunk
% mirror to fill symmetric null matrix
    nEvt = min(size(allUnshift,1), size(allShift,1));
    nAll = min(size(allUnshift,2), size(allShift,2));
    segLen = chunkEnd - chunkStart + 1;

    nullCorrMat = nan(nAll, nAll);

    for i = 1:nAll
        for j = (i+1):nAll
            Avec = nan(1, nEvt * segLen);
            Bvec = nan(1, nEvt * segLen);

            w = 0;
            for e = 1:nEvt
                aseg = squeeze(allUnshift(e, i, chunkStart:chunkEnd));
                bseg = squeeze(allShift(e,   j, chunkStart:chunkEnd));

                idx = (w+1):(w+segLen);
                Avec(idx) = aseg(:)';
                Bvec(idx) = bseg(:)';
                w = w + segLen;
            end

            valid = ~isnan(Avec) & ~isnan(Bvec);
            if nnz(valid) > 10
                r = corr(Avec(valid)', Bvec(valid)');
                nullCorrMat(i,j) = r;
                nullCorrMat(j,i) = r;
            end
        end
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
    ok = false; cortexIntFR = [];
    try
        frFile = fullfile(folderPath, 'NeuralFiringRates1msBins10msGauss.mat');
        load(frFile, 'cortexFRs','cortexInds');

        clsFile = fullfile('/home/asa7288/Transfer', 'AA_classifications.mat');
        load(clsFile, 'classifications');

        if contains(folderPath,'D026'), matchRow=1;
        elseif contains(folderPath,'D020'), matchRow=2;
        elseif contains(folderPath,'D024'), matchRow=3;
        elseif contains(folderPath,'D043'), matchRow=4;
        else, error('unknown session folder'); end

        cls = classifications{matchRow,1}(cortexInds);
        cortexIntFR = cortexFRs(cls==1, :);
        ok = true;
    catch ME
        warning(ME.message);
        ok = false;
    end
end

function [cortexPyrFR, ok] = loadCortexPyrFR(folderPath)
    ok = false; cortexPyrFR = [];
    try
        frFile = fullfile(folderPath, 'NeuralFiringRates1msBins10msGauss.mat');
        load(frFile, 'cortexFRs','cortexInds');

        clsFile = fullfile('/home/asa7288/Transfer', 'AA_classifications.mat');
        load(clsFile, 'classifications');

        if contains(folderPath,'D026'), matchRow=1;
        elseif contains(folderPath,'D020'), matchRow=2;
        elseif contains(folderPath,'D024'), matchRow=3;
        elseif contains(folderPath,'D043'), matchRow=4;
        else, error('unknown session folder'); end

        cls = classifications{matchRow,1}(cortexInds);
        cortexPyrFR = cortexFRs(cls==0, :);
        ok = true;
    catch ME
        warning(ME.message);
        ok = false;
    end
end
