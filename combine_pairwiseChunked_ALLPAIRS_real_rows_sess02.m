% combine_pairwiseChunked_ALLPAIRS_real_rows_sess02.m

% combines per-row real outputs into the standard single file:
%   pairwiseChunkedXCorr_ALLPAIRS_unshifted_sess02_fullxc.mat

% run on quest after all row jobs finish.

clearvars; clc;

%% settings
sessInd = 2;
baseDirs = {
    '/home/asa7288/Transfer/D026', ...
    '/home/asa7288/Transfer/D020', ...
    '/home/asa7288/Transfer/D024', ...
    '/home/asa7288/Transfer/D043'
};
baseDir = baseDirs{sessInd};

outDir = fullfile(baseDir, 'quest_runs');
tmpDir = fullfile(outDir, 'tmp_real_rows');

assert(exist(tmpDir,'dir')==7, 'missing tmp row directory: %s', tmpDir);

% find row files
rowFiles = dir(fullfile(tmpDir, sprintf('pairwiseChunkedXCorr_ALLPAIRS_unshifted_sess%02d_row*.mat', sessInd)));
assert(~isempty(rowFiles), 'no row files found in %s', tmpDir);

% load one file to get dimensions/metadata
one = load(fullfile(tmpDir, rowFiles(1).name));
lags = one.lags;
binSize = one.binSize;
tAxis = one.tAxis;
nAll = one.nAll;

chunkHalf = one.chunkHalf;
channelsToUse = one.channelsToUse;
doBaselineNorm = one.doBaselineNorm;

nPyr_ref = one.nPyr_ref;
nInt_ref = one.nInt_ref;
neuronType = one.neuronType;
pyrIdx = one.pyrIdx;
intIdx = one.intIdx;

nL = numel(lags);

% allocate final matrices (same as original)
xcMat_all = nan(nAll, nAll, nL);
peakCorrMat_all = nan(nAll, nAll);
peakLagSecMat_all = nan(nAll, nAll);

runtimeSec = 0;

% fill rows
fprintf('combining %d row files...\n', numel(rowFiles));
for k = 1:numel(rowFiles)
    f = fullfile(tmpDir, rowFiles(k).name);
    R = load(f, 'jobInd','xcRow','peakCorrRow','peakLagSecRow','runtimeSec');

    i = R.jobInd;

    % only upper-tri (j>i) was computed
    for j = (i+1):nAll
        xcMat_all(i,j,:) = R.xcRow(j,:);
        peakCorrMat_all(i,j) = R.peakCorrRow(j);
        peakLagSecMat_all(i,j) = R.peakLagSecRow(j);
    end

    runtimeSec = runtimeSec + R.runtimeSec;
end

jobInd = 0; % match original "real" job
outFile = fullfile(outDir, sprintf('pairwiseChunkedXCorr_ALLPAIRS_unshifted_sess%02d_fullxc.mat', sessInd));

save(outFile, ...
    'lags','binSize','tAxis', ...
    'xcMat_all','peakCorrMat_all','peakLagSecMat_all', ...
    'chunkHalf','channelsToUse','doBaselineNorm', ...
    'jobInd','baseDir','sessInd', ...
    'nPyr_ref','nInt_ref','nAll','neuronType','pyrIdx','intIdx', ...
    'runtimeSec', ...
    '-v7.3');

fprintf('saved combined file: %s\n', outFile);
