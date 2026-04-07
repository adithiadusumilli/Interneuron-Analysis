% pairwise_chunked_allpairs_REAL_byRow_quest.m
%
% real (unshifted) chunked all-pairs xcorr, parallelized by seed neuron i
% each array task computes one row: i vs j (j>i), full lag sweep
% outputs per-row temp files, later combined into:
%   pairwiseChunkedXCorr_ALLPAIRS_unshifted_sess02_fullxc.mat
%
% required input (from sbatch):
%   jobInd = seed neuron index i (1..nAll)

clc;

%% ---------------- settings ----------------
baseDirs = {
    '/home/asa7288/Transfer/D026', ...
    '/home/asa7288/Transfer/D020', ...
    '/home/asa7288/Transfer/D024', ...
    '/home/asa7288/Transfer/D043'
};

sessInd = 2;
baseDir = baseDirs{sessInd};

channelsToUse   = 1:4;
chunkHalf       = 200;   % 401 samples
maxLagSecs      = 0.5;
doBaselineNorm  = true;

assert(exist('jobInd','var')==1, 'define jobInd before running (1..nAll)');

%% ---------------- load windowed data ----------------
S = load(fullfile(baseDir, 'EMG_Neural_AllChannels.mat'), ...
    'pyrCxWinCell', 'intCxWinCell', 'tAxis');

assert(isfield(S,'tAxis') && isfield(S,'pyrCxWinCell') && isfield(S,'intCxWinCell'), ...
    'missing required vars in EMG_Neural_AllChannels.mat');

tAxis = S.tAxis(:)';     % ms
binSize = 0.001;         % sec
T = numel(tAxis);

[~, cIdx] = min(abs(tAxis - 0));
chunkStart = cIdx - chunkHalf;
chunkEnd   = cIdx + chunkHalf;
assert(chunkStart>=1 && chunkEnd<=T, 'chunkHalf does not fit window length');

reqMaxLagBins = round(maxLagSecs / binSize);
maxLagLeft  = chunkStart - 1;
maxLagRight = T - chunkEnd;
maxLagBins  = min([reqMaxLagBins, maxLagLeft, maxLagRight]);

lags = -maxLagBins:maxLagBins;
nL = numel(lags);

%% ============================================================
% build unshifted pooled windows (exactly like your original)
%% ============================================================
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

    pyrAll_unshift = cat(1, pyrAll_unshift, pyrWin);
    intAll_unshift = cat(1, intAll_unshift, intWin);
end

assert(~isempty(pyrAll_unshift) && ~isempty(intAll_unshift), 'no valid channels found');

nEvt_un = min(size(pyrAll_unshift,1), size(intAll_unshift,1));
pyrAll_unshift = pyrAll_unshift(1:nEvt_un,:,:);
intAll_unshift = intAll_unshift(1:nEvt_un,:,:);

allAll_unshift = cat(2, pyrAll_unshift, intAll_unshift); % events x nAll x time
nAll = size(allAll_unshift,2);

neuronType = [zeros(1,nPyr_ref), ones(1,nInt_ref)]; % 0=pyr, 1=int
pyrIdx = 1:nPyr_ref;
intIdx = (nPyr_ref+1):nAll;

assert(jobInd>=1 && jobInd<=nAll, 'jobInd must be 1..nAll');

fprintf('\n=== REAL ROW JOB | sessInd=%d | jobInd=%d | nEvt=%d | nAll=%d | nLags=%d ===\n', ...
    sessInd, jobInd, nEvt_un, nAll, nL);

%% ---------------- compute only this row ----------------
tic
[xcRow, peakCorrRow, peakLagSecRow] = ...
    lagSweepPairwise_fullxc_chunked_oneRow(allAll_unshift, jobInd, lags, chunkStart, chunkEnd, T, binSize);
runtimeSec = toc;

%% ---------------- save per-row temp file ----------------
outDir = fullfile(baseDir, 'quest_runs');
if ~exist(outDir,'dir'), mkdir(outDir); end

tmpDir = fullfile(outDir, 'tmp_real_rows');
if ~exist(tmpDir,'dir'), mkdir(tmpDir); end

outFile = fullfile(tmpDir, sprintf('pairwiseChunkedXCorr_ALLPAIRS_unshifted_sess%02d_row%03d.mat', sessInd, jobInd));
save(outFile, ...
    'lags','binSize','tAxis', ...
    'xcRow','peakCorrRow','peakLagSecRow', ...
    'chunkHalf','channelsToUse','doBaselineNorm', ...
    'jobInd','baseDir','sessInd', ...
    'nPyr_ref','nInt_ref','nAll','neuronType','pyrIdx','intIdx', ...
    'runtimeSec', ...
    '-v7.3');

fprintf('saved row file: %s (%.1f s)\n', outFile, runtimeSec);

%% ===================== helpers =====================

function Wout = subtractTrialBaseline3d(Win, tAxis, tStart, tEnd)
    [~, iStart] = min(abs(tAxis - tStart));
    [~, iEnd]   = min(abs(tAxis - tEnd));
    if iEnd < iStart, tmp=iStart; iStart=iEnd; iEnd=tmp; end
    base = mean(Win(:,:,iStart:iEnd), 3, 'omitnan');
    Wout = Win - base;
end

function [xcRow, peakCorrRow, peakLagSecRow] = lagSweepPairwise_fullxc_chunked_oneRow(allAll, i, lags, chunkStart, chunkEnd, T, binSize)
    nEvt = size(allAll,1);
    nAll = size(allAll,2);
    nL = numel(lags);
    segLen = chunkEnd - chunkStart + 1;

    xcRow = nan(nAll, nL);
    peakCorrRow = nan(nAll, 1);
    peakLagSecRow = nan(nAll, 1);

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
                xcRow(j,iL) = r;

                if ~isnan(r) && r > bestR
                    bestR = r;
                    bestLagSec = L * binSize;
                end
            end
        end

        if isfinite(bestR)
            peakCorrRow(j) = bestR;
            peakLagSecRow(j) = bestLagSec;
        end
    end
end
