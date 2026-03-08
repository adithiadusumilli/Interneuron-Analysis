function plotNonChunkedXcorrPerBehavior_saved(saveFile)
% plotting-only version for nonchunked behavior xcorr

% j run: plotNonChunkedXcorrPerBehavior_saved('X:\David\AnalysesData\nonchunked_xcorr_by_umap_cortex_allSessions_saved.mat')

arguments
    saveFile (1,1) string
end

S = load(saveFile, 'results');
results = S.results;

nSess = numel(results.sessions);
nB = numel(results.behaviors);

for iDir = 1:nSess
    sess = results.sessions(iDir);

    if isempty(sess.baseDir)
        continue;
    end

    figure('Color','w', 'Name', sprintf('non-chunked xcorr by %s | sess %d', results.labelType, iDir));
    tiledlayout(1, nB, 'TileSpacing','compact','Padding','compact');

    for bIdx = 1:nB
        nexttile; hold on;

        behOut = sess.beh(bIdx);
        beh = behOut.beh;

        if isempty(beh) || behOut.nTimepoints == 0 || isempty(behOut.lagsSec)
            title(sprintf('beh %d (n=0)', results.behaviors(bIdx)));
            axis off;
            continue;
        end

        plot(behOut.lagsSec, behOut.xc, 'k', 'LineWidth', 2);
        xline(0, 'k:');

        if ~any(isnan(behOut.ctrlCorrCI))
            yline(behOut.ctrlCorrCI(1), '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);
            yline(behOut.ctrlCorrCI(2), '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);
        end

        if ~isnan(behOut.peakLagSec)
            xline(behOut.peakLagSec, 'r--', 'LineWidth', 1.5);
        end

        if ~any(isnan(behOut.lagCI))
            xline(behOut.lagCI(1), '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);
            xline(behOut.lagCI(2), '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);
        end

        title(sprintf('beh %d (n=%d)', beh, behOut.nTimepoints));
        xlabel('lag (s)');
        if bIdx == 1
            ylabel('correlation');
        end
        box off;
    end

    sgtitle(sprintf('non-chunked xcorr by %s | maxLag=\\pm%.3f s | %s', results.labelType, results.maxLagSecs, sess.sessionTag));
end
end
