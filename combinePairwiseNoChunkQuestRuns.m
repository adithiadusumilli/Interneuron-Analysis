function combinePairwiseNoChunkQuestRuns()
% combines saved quest outputs from pairwise_nochunk_allPairs_real_or_shift_quest.m

% this function:
%   - loads all real-row files and reconstructs the full symmetric real matrices
%   - loads all 100 shift files and stacks the null correlation matrices
%   - computes 2.5% / 97.5% null percentiles for each pair
%   - saves one combined file per session

% expected save-file naming conventions (from quest code):
%   real:  pairwise_nochunk_allPairs_sess%02d_real_row%03d.mat
%   shift: pairwise_nochunk_allPairs_sess%02d_shift_%03d.mat

% sessions:
%   sess 1 (D026): 54 real rows, 100 shifts
%   sess 2 (D020): 83 real rows, 100 shifts
%   sess 3 (D024): 91 real rows, 100 shifts
%   sess 4 (D043): 46 real rows, 100 shifts

% j run:
%   combinePairwiseNoChunkQuestRuns

clc;

%% ---------- settings ----------
baseDirs = {
    '/home/asa7288/Transfer/D026', ...
    '/home/asa7288/Transfer/D020', ...
    '/home/asa7288/Transfer/D024', ...
    '/home/asa7288/Transfer/D043'
};

nRealRowsExpected = [54, 83, 91, 46];
nShiftsExpected = 100;

for sessInd = 1:numel(baseDirs)
    baseDir = baseDirs{sessInd};
    outDir = fullfile(baseDir, 'quest_runs');

    fprintf('\n==============================\n');
    fprintf('combining session %d: %s\n', sessInd, baseDir);
    fprintf('==============================\n');

    %% ---------- combine real row files ----------
    realPattern = sprintf('pairwise_nochunk_allPairs_sess%02d_real_row*.mat', sessInd);
    realFiles = dir(fullfile(outDir, realPattern));

    if isempty(realFiles)
        warning('no real files found for session %d in %s. skipping session.', sessInd, outDir);
        continue;
    end

    % parse row numbers and sort
    realRows = nan(numel(realFiles),1);
    for k = 1:numel(realFiles)
        tok = regexp(realFiles(k).name, 'real_row(\d+)\.mat$', 'tokens', 'once');
        if ~isempty(tok)
            realRows(k) = str2double(tok{1});
        end
    end
    [realRows, sortIdx] = sort(realRows);
    realFiles = realFiles(sortIdx);

    fprintf('found %d/%d real row files\n', numel(realFiles), nRealRowsExpected(sessInd));

    % load first file to get dimensions
    tmp = load(fullfile(outDir, realFiles(1).name), ...
        'lags', 'binSize', 'nInt', 'nPyr', 'nAll', 'typeVec', 'numBins', 'baseDir');
    lags = tmp.lags;
    binSize = tmp.binSize;
    nInt = tmp.nInt;
    nPyr = tmp.nPyr;
    nAll = tmp.nAll;
    typeVec = tmp.typeVec;
    numBins = tmp.numBins;

    nL = numel(lags);

    % initialize combined real outputs
    xcMatAll = nan(nAll, nAll, nL);
    peakCorrMatAll = nan(nAll, nAll);
    peakLagMatAll = nan(nAll, nAll);
    realRowDone = false(nAll,1);

    for k = 1:numel(realFiles)
        thisFile = fullfile(outDir, realFiles(k).name);
        D = load(thisFile, 'xcRowAll', 'peakCorrRowAll', 'peakLagRowAll', 'jobInd');

        i = D.jobInd;
        realRowDone(i) = true;

        for j = (i+1):nAll
            if j <= size(D.xcRowAll,1) && ~all(isnan(D.xcRowAll(j,:)))
                xcMatAll(i,j,:) = D.xcRowAll(j,:);
                xcMatAll(j,i,:) = D.xcRowAll(j,:);

                peakCorrMatAll(i,j) = D.peakCorrRowAll(j);
                peakCorrMatAll(j,i) = D.peakCorrRowAll(j);

                peakLagMatAll(i,j) = D.peakLagRowAll(j);
                peakLagMatAll(j,i) = D.peakLagRowAll(j);
            end
        end
    end

    % diagonal should remain nan (self-self not computed)
    fprintf('real rows recovered: %d/%d\n', nnz(realRowDone), nAll);

    %% ---------- combine shift files ----------
    shiftPattern = sprintf('pairwise_nochunk_allPairs_sess%02d_shift_*.mat', sessInd);
    shiftFiles = dir(fullfile(outDir, shiftPattern));

    if isempty(shiftFiles)
        warning('no shift files found for session %d in %s. saving real-only combination.', sessInd, outDir);

        combinedFile = fullfile(outDir, sprintf('pairwise_nochunk_allPairs_sess%02d_COMBINED.mat', sessInd));
        save(combinedFile, ...
            'sessInd', 'baseDir', 'lags', 'binSize', 'nInt', 'nPyr', 'nAll', 'typeVec', 'numBins', ...
            'xcMatAll', 'peakCorrMatAll', 'peakLagMatAll', 'realRowDone', ...
            '-v7.3');

        fprintf('saved combined real-only file: %s\n', combinedFile);
        continue;
    end

    % parse shift numbers and sort
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

    %% ---------- compute null percentiles per pair ----------
    nullCorrPrctile2p5 = prctile(nullCorrMatAllShifts, 2.5, 3);
    nullCorrPrctile97p5 = prctile(nullCorrMatAllShifts, 97.5, 3);
    nullCorrMean = mean(nullCorrMatAllShifts, 3, 'omitnan');
    nullCorrStd = std(nullCorrMatAllShifts, 0, 3, 'omitnan');

    %% ---------- save combined session file ----------
    combinedFile = fullfile(outDir, sprintf('pairwise_nochunk_allPairs_sess%02d_COMBINED.mat', sessInd));

    save(combinedFile, ...
        'sessInd', 'baseDir', 'lags', 'binSize', 'nInt', 'nPyr', 'nAll', 'typeVec', 'numBins', ...
        'xcMatAll', 'peakCorrMatAll', 'peakLagMatAll', ...
        'nullCorrMatAllShifts', 'nullCorrPrctile2p5', 'nullCorrPrctile97p5', 'nullCorrMean', 'nullCorrStd', ...
        'realRowDone', 'shiftDone', 'realRows', 'shiftNums', 'shiftAmtPerNeuronAll', ...
        '-v7.3');

    fprintf('saved combined file: %s\n', combinedFile);
end

fprintf('\nall sessions finished.\n');

end
