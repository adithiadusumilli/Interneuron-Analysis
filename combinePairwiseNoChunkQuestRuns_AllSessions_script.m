% combinePairwiseNoChunkQuestRuns_AllSessions_script.m
% combines no-chunk pairwise quest outputs across all 6 sessions
% then saves one master file across sessions

% j run: combinePairwiseNoChunkQuestRuns_AllSessions_script

clc;

%% ---------- settings ----------
baseDirs = {
    '/home/asa7288/Transfer/D026', ...
    '/home/asa7288/Transfer/D020', ...
    '/home/asa7288/Transfer/D024', ...
    '/home/asa7288/Transfer/D043', ...
    '/home/asa7288/Transfer/D050', ...
    '/home/asa7288/Transfer/D054'
};

mouseIDs = {'D026','D020','D024','D043','D050','D054'};

nShiftsExpected = 100;

masterOutFile = '/home/asa7288/pairwise_nochunk_allPairs_ALL_SESSIONS_COMBINED.mat';

%% ---------- initialize master ----------
allSessions = struct();
allSessions.baseDirs = baseDirs;
allSessions.mouseIDs = mouseIDs;

allSessions.sessions = repmat(struct( ...
    'sessInd', [], ...
    'mouseID', '', ...
    'baseDir', '', ...
    'combinedFile', '', ...
    'lags', [], ...
    'binSize', [], ...
    'nInt', [], ...
    'nPyr', [], ...
    'nAll', [], ...
    'typeVec', [], ...
    'numBins', [], ...
    'xcMatAll', [], ...
    'peakCorrMatAll', [], ...
    'peakLagMatAll', [], ...
    'nullCorrMatAllShifts', [], ...
    'nullCorrPrctile2p5', [], ...
    'nullCorrPrctile97p5', [], ...
    'nullCorrMean', [], ...
    'nullCorrStd', [], ...
    'realLoaded', false, ...
    'shiftDone', [], ...
    'shiftNums', [], ...
    'shiftAmtPerNeuronAll', []), 1, numel(baseDirs));

%% ---------- loop through sessions ----------
for sessInd = 1:numel(baseDirs)

    baseDir = baseDirs{sessInd};
    mouseID = mouseIDs{sessInd};
    outDir = fullfile(baseDir, 'quest_runs');

    fprintf('\n==============================\n');
    fprintf('combining session %d: %s\n', sessInd, baseDir);
    fprintf('==============================\n');

    %% ---------- load real full-xcorr file ----------
    realFile = fullfile(outDir, sprintf('pairwise_nochunk_allPairs_sess%02d_real_fullxc.mat', sessInd));

    if ~isfile(realFile)
        warning('missing real file for session %d: %s', sessInd, realFile);
        continue;
    end

    R = load(realFile, ...
        'xcMatAll', 'peakCorrMatAll', 'peakLagMatAll', ...
        'lags', 'binSize', 'nInt', 'nPyr', 'nAll', ...
        'typeVec', 'numBins', 'baseDir');

    xcMatAll = R.xcMatAll;
    peakCorrMatAll = R.peakCorrMatAll;
    peakLagMatAll = R.peakLagMatAll;
    lags = R.lags;
    binSize = R.binSize;
    nInt = R.nInt;
    nPyr = R.nPyr;
    nAll = R.nAll;
    typeVec = R.typeVec;
    numBins = R.numBins;

    fprintf('real file loaded: nInt=%d | nPyr=%d | nAll=%d\n', nInt, nPyr, nAll);

    %% ---------- combine shift files ----------
    shiftPattern = sprintf('pairwise_nochunk_allPairs_sess%02d_shift_*.mat', sessInd);
    shiftFiles = dir(fullfile(outDir, shiftPattern));

    if isempty(shiftFiles)
        warning('no shift files found for session %d in %s', sessInd, outDir);

        nullCorrMatAllShifts = [];
        nullCorrPrctile2p5 = [];
        nullCorrPrctile97p5 = [];
        nullCorrMean = [];
        nullCorrStd = [];
        shiftDone = [];
        shiftNums = [];
        shiftAmtPerNeuronAll = [];
    else
        shiftNums = nan(numel(shiftFiles),1);

        for k = 1:numel(shiftFiles)
            tok = regexp(shiftFiles(k).name, 'shift_(\d+)\.mat$', 'tokens', 'once');
            if ~isempty(tok)
                shiftNums(k) = str2double(tok{1});
            end
        end

        [shiftNums, sortIdx] = sort(shiftNums);
        shiftFiles = shiftFiles(sortIdx);

        fprintf('found %d/%d shift files\n', numel(shiftFiles), nShiftsExpected);

        nullCorrMatAllShifts = nan(nAll, nAll, numel(shiftFiles));
        shiftDone = false(numel(shiftFiles),1);
        shiftAmtPerNeuronAll = nan(nAll, numel(shiftFiles));

        for k = 1:numel(shiftFiles)
            thisFile = fullfile(outDir, shiftFiles(k).name);
            D = load(thisFile, 'nullCorrMatAll', 'jobInd', 'shiftAmtPerNeuron');

            nullCorrMatAllShifts(:,:,k) = D.nullCorrMatAll;
            shiftDone(k) = true;

            if isfield(D, 'shiftAmtPerNeuron') && ~isempty(D.shiftAmtPerNeuron)
                shiftAmtPerNeuronAll(:,k) = D.shiftAmtPerNeuron(:);
            end
        end

        fprintf('shift files recovered: %d/%d\n', nnz(shiftDone), numel(shiftFiles));

        nullCorrPrctile2p5 = prctile(nullCorrMatAllShifts, 2.5, 3);
        nullCorrPrctile97p5 = prctile(nullCorrMatAllShifts, 97.5, 3);
        nullCorrMean = mean(nullCorrMatAllShifts, 3, 'omitnan');
        nullCorrStd = std(nullCorrMatAllShifts, 0, 3, 'omitnan');
    end

    %% ---------- save per-session combined file ----------
    combinedFile = fullfile(outDir, sprintf('pairwise_nochunk_allPairs_sess%02d_COMBINED.mat', sessInd));

    save(combinedFile, ...
        'sessInd', 'mouseID', 'baseDir', ...
        'lags', 'binSize', 'nInt', 'nPyr', 'nAll', 'typeVec', 'numBins', ...
        'xcMatAll', 'peakCorrMatAll', 'peakLagMatAll', ...
        'nullCorrMatAllShifts', 'nullCorrPrctile2p5', 'nullCorrPrctile97p5', ...
        'nullCorrMean', 'nullCorrStd', ...
        'shiftDone', 'shiftNums', 'shiftAmtPerNeuronAll', ...
        '-v7.3');

    fprintf('saved per-session combined file:\n%s\n', combinedFile);

    %% ---------- store in master ----------
    allSessions.sessions(sessInd).sessInd = sessInd;
    allSessions.sessions(sessInd).mouseID = mouseID;
    allSessions.sessions(sessInd).baseDir = baseDir;
    allSessions.sessions(sessInd).combinedFile = combinedFile;

    allSessions.sessions(sessInd).lags = lags;
    allSessions.sessions(sessInd).binSize = binSize;
    allSessions.sessions(sessInd).nInt = nInt;
    allSessions.sessions(sessInd).nPyr = nPyr;
    allSessions.sessions(sessInd).nAll = nAll;
    allSessions.sessions(sessInd).typeVec = typeVec;
    allSessions.sessions(sessInd).numBins = numBins;

    allSessions.sessions(sessInd).xcMatAll = xcMatAll;
    allSessions.sessions(sessInd).peakCorrMatAll = peakCorrMatAll;
    allSessions.sessions(sessInd).peakLagMatAll = peakLagMatAll;

    allSessions.sessions(sessInd).nullCorrMatAllShifts = nullCorrMatAllShifts;
    allSessions.sessions(sessInd).nullCorrPrctile2p5 = nullCorrPrctile2p5;
    allSessions.sessions(sessInd).nullCorrPrctile97p5 = nullCorrPrctile97p5;
    allSessions.sessions(sessInd).nullCorrMean = nullCorrMean;
    allSessions.sessions(sessInd).nullCorrStd = nullCorrStd;

    allSessions.sessions(sessInd).realLoaded = true;
    allSessions.sessions(sessInd).shiftDone = shiftDone;
    allSessions.sessions(sessInd).shiftNums = shiftNums;
    allSessions.sessions(sessInd).shiftAmtPerNeuronAll = shiftAmtPerNeuronAll;
end

%% ---------- summary ----------
fprintf('\n========== master summary ==========\n');

for sessInd = 1:numel(allSessions.sessions)
    if ~allSessions.sessions(sessInd).realLoaded
        fprintf('%s: missing real file\n', mouseIDs{sessInd});
    else
        fprintf('%s: nInt=%d | nPyr=%d | nAll=%d | shifts=%d\n', ...
            mouseIDs{sessInd}, ...
            allSessions.sessions(sessInd).nInt, ...
            allSessions.sessions(sessInd).nPyr, ...
            allSessions.sessions(sessInd).nAll, ...
            numel(allSessions.sessions(sessInd).shiftNums));
    end
end

%% ---------- save master ----------
save(masterOutFile, 'allSessions', '-v7.3');

fprintf('\nsaved master combined file:\n%s\n', masterOutFile);
