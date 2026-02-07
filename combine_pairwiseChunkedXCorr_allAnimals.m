% combine_pairwiseChunkedXCorr_allAnimals.m

% combines pairwise CHUNKED xcorr outputs across animals into one .mat
% - unshifted: full lag sweep (xcMat + peakCorrMat + peakLagSecMat)
% - shifted:   0-lag-only controls saved as nullCorrMat per shift

% saves ONE file to: /home/asa7288/Transfer/

clearvars; clc;

%% ---------------- settings ----------------
baseDirs = {
    '/home/asa7288/Transfer/D026', ...
    '/home/asa7288/Transfer/D020', ...
    '/home/asa7288/Transfer/D024', ...
    '/home/asa7288/Transfer/D043'
};

nShifts = 100;

%% ---------------- holders ----------------
all_xcMat          = cell(numel(baseDirs),1);  % (nInt x nPyr x nLags)
all_peakCorrMat    = cell(numel(baseDirs),1);  % (nInt x nPyr)
all_peakLagSecMat  = cell(numel(baseDirs),1);  % (nInt x nPyr)

all_nullCorrShifts = cell(numel(baseDirs),1);  % (nInt x nPyr x nShifts)

all_lags           = cell(numel(baseDirs),1);
all_binSize        = cell(numel(baseDirs),1);
all_tAxis          = cell(numel(baseDirs),1);

all_chunkHalf      = cell(numel(baseDirs),1);
all_channelsToUse  = cell(numel(baseDirs),1);
all_doBaselineNorm = cell(numel(baseDirs),1);

all_sessInd        = nan(numel(baseDirs),1);   % if present in files
all_baseDirEcho    = cell(numel(baseDirs),1);  % echo what file saved

%% ---------------- main loop ----------------
for iDir = 1:numel(baseDirs)
    baseDir = baseDirs{iDir};
    outDir  = fullfile(baseDir, 'quest_runs');

    fprintf('\n==============================\n');
    fprintf('animal %d: %s\n', iDir, baseDir);

    if ~exist(outDir,'dir')
        warning('missing quest_runs folder: %s (skipping)', outDir);
        continue;
    end

    % ---------------- find unshifted file ----------------
    % support both naming styles:
    %   pairwiseChunkedXCorr_unshifted_fullxc.mat
    %   pairwiseChunkedXCorr_unshifted_sess##_fullxc.mat
    unshiftedCandidates = dir(fullfile(outDir, 'pairwiseChunkedXCorr_unshifted*_fullxc.mat'));

    if isempty(unshiftedCandidates)
        warning('no unshifted fullxc file found in %s (skipping animal)', outDir);
        continue;
    end

    % if multiple match, take the newest
    [~, idxNewest] = max([unshiftedCandidates.datenum]);
    unshiftedFile = fullfile(outDir, unshiftedCandidates(idxNewest).name);
    fprintf('unshifted file: %s\n', unshiftedCandidates(idxNewest).name);

    % ---------------- load unshifted (NO struct) ----------------
    % loads variables directly into workspace
    load(unshiftedFile, ...
        'xcMat','peakCorrMat','peakLagSecMat', ...
        'lags','binSize','tAxis', ...
        'chunkHalf','channelsToUse','doBaselineNorm', ...
        'jobInd','baseDir','sessInd');

    % store unshifted
    all_xcMat{iDir}         = xcMat;
    all_peakCorrMat{iDir}   = peakCorrMat;
    all_peakLagSecMat{iDir} = peakLagSecMat;

    all_lags{iDir}    = lags;
    all_binSize{iDir} = binSize;
    all_tAxis{iDir}   = tAxis;

    if exist('chunkHalf','var'),      all_chunkHalf{iDir}      = chunkHalf;      else, all_chunkHalf{iDir} = []; end
    if exist('channelsToUse','var'),  all_channelsToUse{iDir}  = channelsToUse;  else, all_channelsToUse{iDir} = []; end
    if exist('doBaselineNorm','var'), all_doBaselineNorm{iDir} = doBaselineNorm; else, all_doBaselineNorm{iDir} = []; end

    if exist('sessInd','var') && ~isempty(sessInd)
        all_sessInd(iDir) = sessInd;
    end
    if exist('baseDir','var') && ~isempty(baseDir)
        all_baseDirEcho{iDir} = baseDir;
    else
        all_baseDirEcho{iDir} = baseDirs{iDir};
    end

    % dimensions for shift stacking
    nInt = size(xcMat,1);
    nPyr = size(xcMat,2);

    % clear big vars we no longer need in workspace (optional but keeps things clean)
    clear peakCorrMat peakLagSecMat;

    % ---------------- load shifts and stack ----------------
    nullCorrMat_allShifts = nan(nInt, nPyr, nShifts);

    nFound = 0;
    for s = 1:nShifts
        % support both naming styles:
        %   pairwiseChunkedXCorr_shift_###_zerolag.mat
        %   pairwiseChunkedXCorr_sess##_shift_###_zerolag.mat
        f1 = fullfile(outDir, sprintf('pairwiseChunkedXCorr_shift_%03d_zerolag.mat', s));
        f2 = fullfile(outDir, sprintf('pairwiseChunkedXCorr_sess%02d_shift_%03d_zerolag.mat', iDir, s)); % fallback guess

        if exist(f1,'file')
            shiftFile = f1;
        elseif exist(f2,'file')
            shiftFile = f2;
        else
            % last resort: wildcard search for this shift number
            dd = dir(fullfile(outDir, sprintf('*shift_%03d*_zerolag.mat', s)));
            if isempty(dd)
                warning('missing shift %03d file in %s', s, outDir);
                continue;
            end
            shiftFile = fullfile(outDir, dd(1).name);
        end

        % load nullCorrMat directly
        load(shiftFile, 'nullCorrMat');

        if ~exist('nullCorrMat','var') || isempty(nullCorrMat)
            warning('nullCorrMat missing/empty in %s', shiftFile);
            clear nullCorrMat;
            continue;
        end

        % safety: size check
        if ~isequal(size(nullCorrMat), [nInt, nPyr])
            warning('size mismatch in %s (got %dx%d, expected %dx%d). skipping.', ...
                shiftFile, size(nullCorrMat,1), size(nullCorrMat,2), nInt, nPyr);
            clear nullCorrMat;
            continue;
        end

        nullCorrMat_allShifts(:,:,s) = nullCorrMat;
        nFound = nFound + 1;

        clear nullCorrMat;
    end

    all_nullCorrShifts{iDir} = nullCorrMat_allShifts;

    fprintf('loaded shifts: %d / %d\n', nFound, nShifts);

    % clear remaining unshifted vars to avoid accidental carry-over
    clear xcMat lags binSize tAxis chunkHalf channelsToUse doBaselineNorm jobInd sessInd;
end

%% ---------------- save combined ----------------
outFileAll = '/home/asa7288/Transfer/allAnimals_pairwiseChunkedXCorr_combined.mat';

save(outFileAll, ...
    'all_xcMat','all_peakCorrMat','all_peakLagSecMat', ...
    'all_nullCorrShifts', ...
    'all_lags','all_binSize','all_tAxis', ...
    'all_chunkHalf','all_channelsToUse','all_doBaselineNorm', ...
    'all_sessInd','all_baseDirEcho', ...
    'baseDirs','nShifts', ...
    '-v7.3');

fprintf('\nSAVED combined pairwise chunked results to:\n%s\n', outFileAll);
