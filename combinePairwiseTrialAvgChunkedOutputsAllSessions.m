function combinePairwiseTrialAvgChunkedOutputsAllSessions()

clc;

baseDirs = {
    '/home/asa7288/Transfer/D026', ...
    '/home/asa7288/Transfer/D020', ...
    '/home/asa7288/Transfer/D024', ...
    '/home/asa7288/Transfer/D043'
};

outFile = '/home/asa7288/pairwiseChunkedTrialAvgXCorr_ALLPAIRS_ALLSESS_COMBINED.mat';

nSess = numel(baseDirs);

allSessions = struct();
allSessions.baseDirs = baseDirs;
allSessions.sessions = repmat(struct( ...
    'baseDir', '', ...
    'sessionTag', '', ...
    'lags', [], ...
    'binSize', [], ...
    'tAxis', [], ...
    'chunkHalf', [], ...
    'channelsToUse', [], ...
    'doBaselineNorm', [], ...
    'xcMat_all', [], ...
    'peakCorrMat_all', [], ...
    'peakLagSecMat_all', [], ...
    'nPyr_ref', [], ...
    'nInt_ref', [], ...
    'nAll', [], ...
    'neuronType', [], ...
    'pyrIdx', [], ...
    'intIdx', [], ...
    'realIntPyrLags', [], ...
    'realMedianLag', [], ...
    'nullLagH0', [], ...
    'nullLagH50', [], ...
    'pValH0', [], ...
    'pValH50', [], ...
    'evidenceRatio_H0_over_H50', [], ...
    'nullPeakLagVec_all', [], ...
    'nullPeakCorrVec_all', [], ...
    'nullPairType_all', [], ...
    'perm_inds_found', []), 1, nSess);

for sessInd = 1:nSess

    baseDir = baseDirs{sessInd};
    outDir = fullfile(baseDir, 'quest_runs');
    sessionTag = localShortTag(baseDir);

    fprintf('\n==============================\n');
    fprintf('combining %s | sessInd %d\n', sessionTag, sessInd);
    fprintf('==============================\n');

    sess = struct();
    sess.baseDir = baseDir;
    sess.sessionTag = sessionTag;

    %% ---------- real file ----------
    realFile = fullfile(outDir, ...
        sprintf('pairwiseChunkedTrialAvgXCorr_ALLPAIRS_unshifted_sess%02d_fullxc.mat', sessInd));

    if ~isfile(realFile)
        warning('missing real file: %s', realFile);
    else
        R = load(realFile);

        sess.lags = R.lags;
        sess.binSize = R.binSize;
        sess.tAxis = R.tAxis;
        sess.chunkHalf = R.chunkHalf;
        sess.channelsToUse = R.channelsToUse;
        sess.doBaselineNorm = R.doBaselineNorm;

        sess.xcMat_all = R.xcMat_all;
        sess.peakCorrMat_all = R.peakCorrMat_all;
        sess.peakLagSecMat_all = R.peakLagSecMat_all;

        sess.nPyr_ref = R.nPyr_ref;
        sess.nInt_ref = R.nInt_ref;
        sess.nAll = R.nAll;
        sess.neuronType = R.neuronType;
        sess.pyrIdx = R.pyrIdx;
        sess.intIdx = R.intIdx;

        % Bayes outputs, if already appended
        bayesFields = { ...
            'realIntPyrLags', ...
            'realMedianLag', ...
            'nullLagH0', ...
            'nullLagH50', ...
            'pValH0', ...
            'pValH50', ...
            'evidenceRatio_H0_over_H50'};

        for f = 1:numel(bayesFields)
            thisField = bayesFields{f};

            if isfield(R, thisField)
                sess.(thisField) = R.(thisField);
            else
                sess.(thisField) = [];
            end
        end
    end

    %% ---------- null files ----------
    nullFiles = dir(fullfile(outDir, ...
        sprintf('pairwiseChunkedTrialAvgXCorr_ALLPAIRS_sess%02d_permNull_*.mat', sessInd)));

    nullPeakLagVec_all = [];
    nullPeakCorrVec_all = [];
    nullPairType_all = strings(0,1);
    perm_inds_found = [];

    if isempty(nullFiles)
        warning('no null files found for %s', baseDir);
    else
        nullNums = nan(numel(nullFiles),1);

        for k = 1:numel(nullFiles)
            tok = regexp(nullFiles(k).name, 'permNull_(\d+)\.mat$', 'tokens', 'once');
            if ~isempty(tok)
                nullNums(k) = str2double(tok{1});
            end
        end

        [nullNums, sortIdx] = sort(nullNums);
        nullFiles = nullFiles(sortIdx);

        for k = 1:numel(nullFiles)
            D = load(fullfile(outDir, nullFiles(k).name), ...
                'nullPeakLagVec', 'nullPeakCorrVec', 'nullPairType');

            if isfield(D, 'nullPeakLagVec')
                nullPeakLagVec_all = [nullPeakLagVec_all; D.nullPeakLagVec(:)];
            end

            if isfield(D, 'nullPeakCorrVec')
                nullPeakCorrVec_all = [nullPeakCorrVec_all; D.nullPeakCorrVec(:)];
            end

            if isfield(D, 'nullPairType')
                nullPairType_all = [nullPairType_all; string(D.nullPairType(:))];
            end
        end

        perm_inds_found = nullNums;
    end

    sess.nullPeakLagVec_all = nullPeakLagVec_all;
    sess.nullPeakCorrVec_all = nullPeakCorrVec_all;
    sess.nullPairType_all = nullPairType_all;
    sess.perm_inds_found = perm_inds_found;

    allSessions.sessions(sessInd) = sess;

    fprintf('real loaded: %d\n', ~isempty(sess.peakLagSecMat_all));
    fprintf('null jobs found: %d\n', numel(sess.perm_inds_found));

    if ~isempty(sess.evidenceRatio_H0_over_H50)
        fprintf('BF H0/H50 loaded: %.6f\n', sess.evidenceRatio_H0_over_H50);
    else
        fprintf('BF H0/H50 not found\n');
    end
end

save(outFile, 'allSessions', '-v7.3');
fprintf('\nsaved combined file:\n%s\n', outFile);

end

function tag = localShortTag(baseDir)
    [~, folderName] = fileparts(baseDir);
    m = regexp(folderName, 'D\d+', 'match', 'once');

    if isempty(m)
        tag = folderName;
    else
        tag = m;
    end
end
