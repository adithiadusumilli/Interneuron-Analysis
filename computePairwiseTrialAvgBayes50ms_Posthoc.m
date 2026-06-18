function computePairwiseTrialAvgBayes50ms_Posthoc()
% pairwise bayes-style evidence ratio
% H0: real int-pyr median peak lag compared to all-vs-all random null
% H50: same null shifted by +0.050 s
% pos lag = interneuron leads pyramidal

baseDirs = {
    '/home/asa7288/Transfer/D026', ...
    '/home/asa7288/Transfer/D020', ...
    '/home/asa7288/Transfer/D024', ...
    '/home/asa7288/Transfer/D043'
};

outSummary = struct();

for sessInd = 1:numel(baseDirs)

    baseDir = baseDirs{sessInd};
    outDir = fullfile(baseDir, 'quest_runs');

    realFile = fullfile(outDir, ...
        sprintf('pairwiseChunkedTrialAvgXCorr_ALLPAIRS_unshifted_sess%02d_fullxc.mat', sessInd));

    if ~isfile(realFile)
        warning('missing real file: %s', realFile);
        continue;
    end

    R = load(realFile, 'peakLagSecMat_all','peakCorrMat_all','pyrIdx','intIdx','baseDir');

    realIntPyrLags = R.peakLagSecMat_all(R.pyrIdx, R.intIdx);
    realIntPyrLags = realIntPyrLags(:);
    realIntPyrLags = realIntPyrLags(~isnan(realIntPyrLags) & isfinite(realIntPyrLags));

    realMedianLag = median(realIntPyrLags, 'omitnan');

    nullFiles = dir(fullfile(outDir, sprintf('pairwiseChunkedTrialAvgXCorr_ALLPAIRS_sess%02d_permNull_*.mat', sessInd)));

    nullLagH0 = [];

    for k = 1:numel(nullFiles)
        D = load(fullfile(outDir, nullFiles(k).name), 'nullPeakLagVec');
        nullLagH0 = [nullLagH0; D.nullPeakLagVec(:)];
    end

    nullLagH0 = nullLagH0(~isnan(nullLagH0) & isfinite(nullLagH0));
    nullLagH50 = nullLagH0 + 0.050;

    realDistH0 = abs(realMedianLag - 0);
    nullDistH0 = abs(nullLagH0 - 0);

    realDistH50 = abs(realMedianLag - 0.050);
    nullDistH50 = abs(nullLagH50 - 0.050);

    pValH0 = (sum(nullDistH0 >= realDistH0) + 1) / (numel(nullDistH0) + 1);
    pValH50 = (sum(nullDistH50 >= realDistH50) + 1) / (numel(nullDistH50) + 1);

    evidenceRatio_H0_over_H50 = pValH0 / pValH50;

    save(realFile, ...
        'realIntPyrLags','realMedianLag', ...
        'nullLagH0','nullLagH50', ...
        'pValH0','pValH50','evidenceRatio_H0_over_H50', ...
        '-append');

    outSummary(sessInd).baseDir = baseDir;
    outSummary(sessInd).realMedianLag = realMedianLag;
    outSummary(sessInd).pValH0 = pValH0;
    outSummary(sessInd).pValH50 = pValH50;
    outSummary(sessInd).evidenceRatio_H0_over_H50 = evidenceRatio_H0_over_H50;

    fprintf('\n%s\n', baseDir);
    fprintf('real median int-pyr lag = %.6f s\n', realMedianLag);
    fprintf('pValH0 = %.6f | pValH50 = %.6f | ratio H0/H50 = %.6f\n', ...
        pValH0, pValH50, evidenceRatio_H0_over_H50);
end

save('/home/asa7288/pairwiseTrialAvgBayes50ms_summary.mat', 'outSummary', '-v7.3');

end
