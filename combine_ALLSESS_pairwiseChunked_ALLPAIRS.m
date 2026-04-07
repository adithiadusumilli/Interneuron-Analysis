% combine_ALLSESS_pairwiseChunked_ALLPAIRS.m

% combines, for each baseDir/session:
%   - the final REAL file:
%       pairwiseChunkedXCorr_ALLPAIRS_unshifted_sess##_fullxc.mat 
%       -- this is after using the combine_pairwiseChunked_ALLPAIRS_real_rows_sess02.m script
%   - the 100 SHIFT files:
%       pairwiseChunkedXCorr_ALLPAIRS_sess##_shift_###_zerolag.mat
%
% into ONE combined file across all sessions.
%
% this assumes sess 2 real has already been reconstructed into the same
% final naming convention as the other sessions.

clearvars; clc;

%% ---------------- settings ----------------
baseDirs = {
    '/home/asa7288/Transfer/D026', ...
    '/home/asa7288/Transfer/D020', ...
    '/home/asa7288/Transfer/D024', ...
    '/home/asa7288/Transfer/D043'
};

sessNames = {'D026','D020','D024','D043'};
nShifts = 100;

outRoot = '/home/asa7288/Transfer';
outFile = fullfile(outRoot, 'pairwiseChunked_ALLPAIRS_ALLSESS_REAL_AND_SHIFTS_COMBINED.mat');

%% ---------------- holders ----------------
all_xcMat_all = cell(numel(baseDirs),1);          % each: nAll x nAll x nL
all_peakCorrMat_all = cell(numel(baseDirs),1);    % each: nAll x nAll
all_peakLagSecMat_all = cell(numel(baseDirs),1);  % each: nAll x nAll

all_nullCorrMat_allShifts = cell(numel(baseDirs),1); % each: nAll x nAll x nShifts

all_lags = cell(numel(baseDirs),1);
all_binSize = cell(numel(baseDirs),1);
all_tAxis = cell(numel(baseDirs),1);

all_chunkHalf = cell(numel(baseDirs),1);
all_channelsToUse = cell(numel(baseDirs),1);
all_doBaselineNorm = cell(numel(baseDirs),1);

all_nPyr_ref = nan(numel(baseDirs),1);
all_nInt_ref = nan(numel(baseDirs),1);
all_nAll = nan(numel(baseDirs),1);

all_neuronType = cell(numel(baseDirs),1);
all_pyrIdx = cell(numel(baseDirs),1);
all_intIdx = cell(numel(baseDirs),1);

all_runtimeRealSec = nan(numel(baseDirs),1);
all_runtimeShiftSec = cell(numel(baseDirs),1);   % each: 1 x nShifts

realFound = false(numel(baseDirs),1);
shiftFound = false(numel(baseDirs), nShifts);

%% ---------------- main loop ----------------
for s = 1:numel(baseDirs)
    baseDir = baseDirs{s};
    outDir = fullfile(baseDir, 'quest_runs');

    fprintf('\n==============================\n');
    fprintf('session %d: %s\n', s, baseDir);

    if ~exist(outDir,'dir')
        warning('missing quest_runs folder: %s', outDir);
        continue;
    end

    %% ---------- load REAL ----------
    realFile = fullfile(outDir, sprintf('pairwiseChunkedXCorr_ALLPAIRS_unshifted_sess%02d_fullxc.mat', s));

    if ~exist(realFile, 'file')
        warning('missing real file: %s', realFile);
        continue;
    end

    R = load(realFile, ...
        'lags','binSize','tAxis', ...
        'xcMat_all','peakCorrMat_all','peakLagSecMat_all', ...
        'chunkHalf','channelsToUse','doBaselineNorm', ...
        'jobInd','baseDir','sessInd', ...
        'nPyr_ref','nInt_ref','nAll','neuronType','pyrIdx','intIdx', ...
        'runtimeSec');

    all_xcMat_all{s} = R.xcMat_all;
    all_peakCorrMat_all{s} = R.peakCorrMat_all;
    all_peakLagSecMat_all{s} = R.peakLagSecMat_all;

    all_lags{s} = R.lags;
    all_binSize{s} = R.binSize;
    all_tAxis{s} = R.tAxis;

    all_chunkHalf{s} = R.chunkHalf;
    all_channelsToUse{s} = R.channelsToUse;
    all_doBaselineNorm{s} = R.doBaselineNorm;

    all_nPyr_ref(s) = R.nPyr_ref;
    all_nInt_ref(s) = R.nInt_ref;
    all_nAll(s) = R.nAll;

    all_neuronType{s} = R.neuronType;
    all_pyrIdx{s} = R.pyrIdx;
    all_intIdx{s} = R.intIdx;

    if isfield(R,'runtimeSec')
        all_runtimeRealSec(s) = R.runtimeSec;
    end

    realFound(s) = true;

    fprintf('loaded REAL: nAll=%d, nLags=%d\n', R.nAll, numel(R.lags));

    %% ---------- load SHIFTS ----------
    nAll = R.nAll;
    nullCorrMat_allShifts = nan(nAll, nAll, nShifts);
    runtimeShiftVec = nan(1, nShifts);

    for sh = 1:nShifts
        shiftFile = fullfile(outDir, sprintf('pairwiseChunkedXCorr_ALLPAIRS_sess%02d_shift_%03d_zerolag.mat', s, sh));

        if ~exist(shiftFile, 'file')
            warning('missing shift file: %s', shiftFile);
            continue;
        end

        S = load(shiftFile, 'nullCorrMat_all', 'runtimeSec');

        if ~isfield(S,'nullCorrMat_all') || isempty(S.nullCorrMat_all)
            warning('nullCorrMat_all missing in %s', shiftFile);
            continue;
        end

        if ~isequal(size(S.nullCorrMat_all), [nAll, nAll])
            warning('size mismatch in %s, got %s expected [%d %d]', ...
                shiftFile, mat2str(size(S.nullCorrMat_all)), nAll, nAll);
            continue;
        end

        nullCorrMat_allShifts(:,:,sh) = S.nullCorrMat_all;
        shiftFound(s,sh) = true;

        if isfield(S,'runtimeSec')
            runtimeShiftVec(sh) = S.runtimeSec;
        end
    end

    all_nullCorrMat_allShifts{s} = nullCorrMat_allShifts;
    all_runtimeShiftSec{s} = runtimeShiftVec;

    fprintf('loaded shifts: %d / %d\n', nnz(shiftFound(s,:)), nShifts);
end

%% ---------------- save combined ----------------
save(outFile, ...
    'baseDirs','sessNames','nShifts', ...
    'all_xcMat_all','all_peakCorrMat_all','all_peakLagSecMat_all', ...
    'all_nullCorrMat_allShifts', ...
    'all_lags','all_binSize','all_tAxis', ...
    'all_chunkHalf','all_channelsToUse','all_doBaselineNorm', ...
    'all_nPyr_ref','all_nInt_ref','all_nAll', ...
    'all_neuronType','all_pyrIdx','all_intIdx', ...
    'all_runtimeRealSec','all_runtimeShiftSec', ...
    'realFound','shiftFound', ...
    '-v7.3');

fprintf('\nSAVED combined file:\n%s\n', outFile);

%% ---------------- summary ----------------
fprintf('\nsummary:\n');
for s = 1:numel(baseDirs)
    fprintf('sess %02d (%s): real=%d, shifts=%d/%d\n', ...
        s, sessNames{s}, realFound(s), nnz(shiftFound(s,:)), nShifts);
end
