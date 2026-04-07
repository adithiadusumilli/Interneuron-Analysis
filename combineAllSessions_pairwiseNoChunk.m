function combineAllSessions_pairwiseNoChunk()
% combines the per-session COMBINED files from combinePairwiseNoChunkQuestRuns.m into one master file
% and saves it into /home/asa7288

% expected inputs (already created by combinePairwiseNoChunkQuestRuns):
%   /home/asa7288/Transfer/D026/quest_runs/pairwise_nochunk_allPairs_sess01_COMBINED.mat
%   /home/asa7288/Transfer/D020/quest_runs/pairwise_nochunk_allPairs_sess02_COMBINED.mat
%   /home/asa7288/Transfer/D024/quest_runs/pairwise_nochunk_allPairs_sess03_COMBINED.mat
%   /home/asa7288/Transfer/D043/quest_runs/pairwise_nochunk_allPairs_sess04_COMBINED.mat

% output:
%   /home/asa7288/pairwise_nochunk_allPairs_ALL_SESSIONS_COMBINED.mat

% j run:
%   combineAllSessions_pairwiseNoChunk

clc;

%% ---------- settings ----------
baseDirs = {
    '/home/asa7288/Transfer/D026', ...
    '/home/asa7288/Transfer/D020', ...
    '/home/asa7288/Transfer/D024', ...
    '/home/asa7288/Transfer/D043'
};

outRoot = '/home/asa7288';
masterOutFile = fullfile(outRoot, 'pairwise_nochunk_allPairs_ALL_SESSIONS_COMBINED.mat');

%% ---------- initialize master struct ----------
allSessions = struct();
allSessions.baseDirs = baseDirs;
allSessions.sessions = repmat(struct( ...
    'sessInd', [], ...
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
    'realRowDone', [], ...
    'shiftDone', [], ...
    'realRows', [], ...
    'shiftNums', [], ...
    'shiftAmtPerNeuronAll', []), 1, numel(baseDirs));

%% ---------- load each session combined file ----------
for sessInd = 1:numel(baseDirs)
    baseDir = baseDirs{sessInd};
    inFile = fullfile(baseDir, 'quest_runs', sprintf('pairwise_nochunk_allPairs_sess%02d_COMBINED.mat', sessInd));

    fprintf('\nloading session %d: %s\n', sessInd, inFile);

    if ~isfile(inFile)
        warning('missing combined file for session %d: %s', sessInd, inFile);
        continue;
    end

    D = load(inFile);

    allSessions.sessions(sessInd).sessInd = sessInd;
    allSessions.sessions(sessInd).baseDir = baseDir;
    allSessions.sessions(sessInd).combinedFile = inFile;

    % required / common fields
    if isfield(D, 'lags'), allSessions.sessions(sessInd).lags = D.lags; end
    if isfield(D, 'binSize'), allSessions.sessions(sessInd).binSize = D.binSize; end
    if isfield(D, 'nInt'), allSessions.sessions(sessInd).nInt = D.nInt; end
    if isfield(D, 'nPyr'), allSessions.sessions(sessInd).nPyr = D.nPyr; end
    if isfield(D, 'nAll'), allSessions.sessions(sessInd).nAll = D.nAll; end
    if isfield(D, 'typeVec'), allSessions.sessions(sessInd).typeVec = D.typeVec; end
    if isfield(D, 'numBins'), allSessions.sessions(sessInd).numBins = D.numBins; end

    % real outputs
    if isfield(D, 'xcMatAll'), allSessions.sessions(sessInd).xcMatAll = D.xcMatAll; end
    if isfield(D, 'peakCorrMatAll'), allSessions.sessions(sessInd).peakCorrMatAll = D.peakCorrMatAll; end
    if isfield(D, 'peakLagMatAll'), allSessions.sessions(sessInd).peakLagMatAll = D.peakLagMatAll; end

    % null outputs
    if isfield(D, 'nullCorrMatAllShifts'), allSessions.sessions(sessInd).nullCorrMatAllShifts = D.nullCorrMatAllShifts; end
    if isfield(D, 'nullCorrPrctile2p5'), allSessions.sessions(sessInd).nullCorrPrctile2p5 = D.nullCorrPrctile2p5; end
    if isfield(D, 'nullCorrPrctile97p5'), allSessions.sessions(sessInd).nullCorrPrctile97p5 = D.nullCorrPrctile97p5; end
    if isfield(D, 'nullCorrMean'), allSessions.sessions(sessInd).nullCorrMean = D.nullCorrMean; end
    if isfield(D, 'nullCorrStd'), allSessions.sessions(sessInd).nullCorrStd = D.nullCorrStd; end

    % bookkeeping
    if isfield(D, 'realRowDone'), allSessions.sessions(sessInd).realRowDone = D.realRowDone; end
    if isfield(D, 'shiftDone'), allSessions.sessions(sessInd).shiftDone = D.shiftDone; end
    if isfield(D, 'realRows'), allSessions.sessions(sessInd).realRows = D.realRows; end
    if isfield(D, 'shiftNums'), allSessions.sessions(sessInd).shiftNums = D.shiftNums; end
    if isfield(D, 'shiftAmtPerNeuronAll'), allSessions.sessions(sessInd).shiftAmtPerNeuronAll = D.shiftAmtPerNeuronAll; end

    fprintf('  loaded: nInt=%d | nPyr=%d | nAll=%d\n', ...
        allSessions.sessions(sessInd).nInt, ...
        allSessions.sessions(sessInd).nPyr, ...
        allSessions.sessions(sessInd).nAll);
end

%% ---------- optional summary ----------
fprintf('\n========== master summary ==========\n');
for sessInd = 1:numel(allSessions.sessions)
    if isempty(allSessions.sessions(sessInd).baseDir)
        fprintf('session %d: missing\n', sessInd);
    else
        fprintf('session %d: nInt=%d | nPyr=%d | nAll=%d\n', ...
            sessInd, ...
            allSessions.sessions(sessInd).nInt, ...
            allSessions.sessions(sessInd).nPyr, ...
            allSessions.sessions(sessInd).nAll);
    end
end

%% ---------- save master file ----------
save(masterOutFile, 'allSessions', '-v7.3');
fprintf('\nsaved master combined file:\n%s\n', masterOutFile);

end
