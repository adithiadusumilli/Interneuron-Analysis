function combineTrialAvgChunkedPopavgOutputsAllSessions()
% combines the trial-averaged chunked population-averaged cross-correlation output files
% across all 6 sessions into one master .mat file

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

mouseIDs = {'D026','D020','D024','D043','D050','D054'};

outFile = '/home/asa7288/concatCrossCorr_trialavg_chunked_popavg_ALL_SESSIONS.mat';

nSess = numel(baseDirs);

%% ---------------- initialize master struct ----------------
allSessions = struct();
allSessions.baseDirs = baseDirs;
allSessions.mouseIDs = mouseIDs;

allSessions.sessions = repmat(struct( ...
    'baseDir', '', ...
    'mouseID', '', ...
    'sessionTag', '', ...
    'lags', [], ...
    'binSize', [], ...
    'chunkHalf', [], ...
    'channelsToUse', [], ...
    'channelsToUseThisSess', [], ...
    'doBaselineNorm', [], ...
    'real_xc', [], ...
    'real_peakLagSec', [], ...
    'real_pyrAvgTrace', [], ...
    'real_intAvgTrace', [], ...
    'real_nTrialsUsed', [], ...
    'perm_xc', [], ...
    'perm_peakLagSec', [], ...
    'perm_nTrialsUsed', [], ...
    'perm_inds_found', [], ...
    'shift_xcZeroLag', [], ...
    'shift_nTrialsUsed', [], ...
    'shift_inds_found', []), 1, nSess);

%% ---------------- loop over sessions ----------------
for iDir = 1:nSess
    baseDir = baseDirs{iDir};
    outDir = fullfile(baseDir, 'quest_runs');

    fprintf('\n==============================\n');
    fprintf('combining session %d: %s\n', iDir, baseDir);
    fprintf('==============================\n');

    sess = allSessions.sessions(iDir);
    sess.baseDir = baseDir;
    sess.mouseID = mouseIDs{iDir};
    sess.sessionTag = localShortTag(baseDir);

    %% ---------- load unpermuted real ----------
    realFile = fullfile(outDir, 'concatCrossCorr_trialavg_chunked_popavg_unperm.mat');

    if ~isfile(realFile)
        warning('missing real file for %s', baseDir);
    else
        R = load(realFile, ...
            'lags','binSize','xc','peakLagSec','chunkHalf','channelsToUse','channelsToUseThisSess', ...
            'doBaselineNorm','pyrAvgTrace','intAvgTrace','nTrialsUsed');

        sess.lags = R.lags;
        sess.binSize = R.binSize;
        sess.chunkHalf = R.chunkHalf;
        sess.channelsToUse = getFieldOrEmpty(R,'channelsToUse');
        sess.channelsToUseThisSess = getFieldOrEmpty(R,'channelsToUseThisSess');
        sess.doBaselineNorm = R.doBaselineNorm;
        sess.real_xc = R.xc;
        sess.real_peakLagSec = R.peakLagSec;
        sess.real_pyrAvgTrace = R.pyrAvgTrace;
        sess.real_intAvgTrace = R.intAvgTrace;
        sess.real_nTrialsUsed = R.nTrialsUsed;
    end

    %% ---------- load permutation files ----------
    permFiles = dir(fullfile(outDir, 'concatCrossCorr_trialavg_chunked_popavg_perm_*.mat'));

    if isempty(permFiles)
        warning('no permutation files found for %s', baseDir);
        sess.perm_xc = [];
        sess.perm_peakLagSec = [];
        sess.perm_nTrialsUsed = [];
        sess.perm_inds_found = [];
    else
        permNums = nan(numel(permFiles),1);

        for k = 1:numel(permFiles)
            tok = regexp(permFiles(k).name, 'perm_(\d+)\.mat$', 'tokens', 'once');
            if ~isempty(tok)
                permNums(k) = str2double(tok{1});
            end
        end

        [permNums, sortIdx] = sort(permNums);
        permFiles = permFiles(sortIdx);

        if isempty(sess.lags)
            tmp = load(fullfile(outDir, permFiles(1).name), ...
                'lags','binSize','chunkHalf','channelsToUse','channelsToUseThisSess','doBaselineNorm');

            sess.lags = tmp.lags;
            sess.binSize = tmp.binSize;
            sess.chunkHalf = tmp.chunkHalf;
            sess.channelsToUse = getFieldOrEmpty(tmp,'channelsToUse');
            sess.channelsToUseThisSess = getFieldOrEmpty(tmp,'channelsToUseThisSess');
            sess.doBaselineNorm = tmp.doBaselineNorm;
        end

        nPerm = numel(permFiles);
        nLags = numel(sess.lags);

        perm_xc = nan(nPerm, nLags);
        perm_peakLagSec = nan(nPerm, 1);
        perm_nTrialsUsed = nan(nPerm, 1);

        for k = 1:nPerm
            D = load(fullfile(outDir, permFiles(k).name), ...
                'xc', 'peakLagSec', 'nTrialsUsed');

            perm_xc(k,:) = D.xc;
            perm_peakLagSec(k) = D.peakLagSec;

            if isfield(D, 'nTrialsUsed')
                perm_nTrialsUsed(k) = D.nTrialsUsed;
            end
        end

        sess.perm_xc = perm_xc;
        sess.perm_peakLagSec = perm_peakLagSec;
        sess.perm_nTrialsUsed = perm_nTrialsUsed;
        sess.perm_inds_found = permNums;
    end

    %% ---------- load shifted zero-lag files ----------
    shiftFiles = dir(fullfile(outDir, 'concatCrossCorr_trialavg_chunked_popavg_shift_*_zerolag.mat'));

    if isempty(shiftFiles)
        warning('no shifted-control files found for %s', baseDir);
        sess.shift_xcZeroLag = [];
        sess.shift_nTrialsUsed = [];
        sess.shift_inds_found = [];
    else
        shiftNums = nan(numel(shiftFiles),1);

        for k = 1:numel(shiftFiles)
            tok = regexp(shiftFiles(k).name, 'shift_(\d+)_zerolag\.mat$', 'tokens', 'once');
            if ~isempty(tok)
                shiftNums(k) = str2double(tok{1});
            end
        end

        [shiftNums, sortIdx] = sort(shiftNums);
        shiftFiles = shiftFiles(sortIdx);

        nShift = numel(shiftFiles);
        shift_xcZeroLag = nan(nShift, 1);
        shift_nTrialsUsed = nan(nShift, 1);

        for k = 1:nShift
            D = load(fullfile(outDir, shiftFiles(k).name), ...
                'xcZeroLag', 'nTrialsUsed');

            shift_xcZeroLag(k) = D.xcZeroLag;

            if isfield(D, 'nTrialsUsed')
                shift_nTrialsUsed(k) = D.nTrialsUsed;
            end
        end

        sess.shift_xcZeroLag = shift_xcZeroLag;
        sess.shift_nTrialsUsed = shift_nTrialsUsed;
        sess.shift_inds_found = shiftNums;
    end

    %% ---------- store ----------
    allSessions.sessions(iDir) = sess;

    fprintf('animal: %s\n', sess.mouseID);
    fprintf('real loaded: %d\n', ~isempty(sess.real_xc));
    fprintf('perm files found: %d\n', numel(sess.perm_inds_found));
    fprintf('shift files found: %d\n', numel(sess.shift_inds_found));
end

%% ---------------- save ----------------
save(outFile, 'allSessions', '-v7.3');
fprintf('\nsaved combined master file:\n%s\n', outFile);

end

%% ================= helpers =================
function tag = localShortTag(baseDir)
    [~, folderName] = fileparts(baseDir);
    m = regexp(folderName, 'D\d+', 'match', 'once');

    if isempty(m)
        tag = folderName;
    else
        tag = m;
    end
end

function val = getFieldOrEmpty(S, fieldName)
    if isfield(S, fieldName)
        val = S.(fieldName);
    else
        val = [];
    end
end
