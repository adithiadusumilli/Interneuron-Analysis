function plotSavedNoChunkXCByBehavior_perAnimal(noChunkFile, minBouts)
% plots saved nonchunked popavg xc results
% split by classifier behv

% for each animal, makes 2 figures:
%   1. lag vs correlation, with 1 panel per behv that passed minBouts
%   2. permutation lag histograms, with 1 panel per behv that passed minBouts
%
% reminder!!!! --> a "bout" / "epoch" is one contiguous stretch of timepoints in timeIdx

% j run: plotSavedNoChunkXCByBehavior_perAnimal("X:\David\AnalysesData\nonchunked_xcorr_by_classifier_cortex_allSessions_saved.mat", 100)

arguments
    noChunkFile (1,1) string
    minBouts (1,1) double = 100
end

behNames = {'unlabeled','climbdown','climbup','eating','grooming','jumpdown','jumping','rearing','still','walkflat','walkgrid'};

origColor = [0 0 0];
corrCIColor = [0 0.2 0.6];
peakLagColor = [0.95 0.45 0.35];
lagCIColor = [0 0 0];

A = load(noChunkFile, 'results');
R = A.results;

nSess = numel(R.sessions);
animalIDs = cell(1, nSess);

for i = 1:nSess
    a = regexp(R.sessions(i).baseDir, 'D\d+', 'match', 'once');
    if isempty(a)
        a = sprintf('Session %d', i);
    end
    animalIDs{i} = a;
end

binSize = 0.001; % nonchunked file was generated at 1 ms bins

fprintf('\n================ nonchunked per animal ================\n')

for s = 1:nSess
    goodBehLabels = {};
    goodBehStruct = {};
    goodPermLags = {};
    goodNBouts = [];
    goodDurSec = [];

    for bIdx = 1:numel(R.behaviors)
        beh = R.behaviors(bIdx);

        if isempty(R.sessions(s).beh)
            continue
        end

        thisBeh = R.sessions(s).beh(bIdx);

        if isempty(thisBeh) || ~isfield(thisBeh, 'timeIdx') || isempty(thisBeh.timeIdx)
            continue
        end

        if isempty(thisBeh.xc) || isempty(thisBeh.lagsSec)
            continue
        end

        timeIdx = sort(thisBeh.timeIdx(:));

        if isempty(timeIdx)
            nBouts = 0;
        else
            nBouts = 1 + sum(diff(timeIdx) > 1);
        end

        if nBouts < minBouts
            continue
        end

        durSec = numel(timeIdx) * binSize;

        if beh >= 0 && beh <= 10
            goodBehLabels{end+1} = sprintf('%d: %s', beh, behNames{beh+1}); %#ok<AGROW>
        else
            goodBehLabels{end+1} = sprintf('%d', beh); %#ok<AGROW>
        end

        goodBehStruct{end+1} = thisBeh; %#ok<AGROW>
        goodPermLags{end+1} = thisBeh.permPeakLags; %#ok<AGROW>
        goodNBouts(end+1) = nBouts; %#ok<AGROW>
        goodDurSec(end+1) = durSec; %#ok<AGROW>
    end

    fprintf('%s nonchunk: %d behaviors passed minBouts=%d\n', ...
        animalIDs{s}, numel(goodBehLabels), minBouts)

    if isempty(goodBehLabels)
        continue
    end

    nPlot = numel(goodBehLabels);
    nCols = min(4, nPlot);
    nRows = ceil(nPlot / nCols);

    % ---- figure 1: lag vs correlation ----
    figure('Name', sprintf('%s nonchunk lag vs correlation by behavior', animalIDs{s}), 'Color', 'w');
    tile_lay = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

    for k = 1:nPlot
        nexttile(tile_lay, k); hold on;

        thisBeh = goodBehStruct{k};
        lagsSec = thisBeh.lagsSec(:)';
        xc = thisBeh.xc(:)';
        peakLag = thisBeh.peakLagSec;
        corrCI = thisBeh.ctrlCorrCI;
        lagCI = thisBeh.lagCI;

        hOrig = plot(lagsSec, xc, 'Color', origColor, 'LineWidth', 2);

        hCorr25 = gobjects(1);
        hCorr97 = gobjects(1);
        if ~isempty(corrCI) && ~any(isnan(corrCI))
            hCorr25 = yline(corrCI(1), '--', 'Color', corrCIColor, 'LineWidth', 1.4);
            hCorr97 = yline(corrCI(2), '--', 'Color', corrCIColor, 'LineWidth', 1.4);
        end

        hPeak = xline(peakLag, '--', 'Color', peakLagColor, 'LineWidth', 1.8);

        hLag1 = gobjects(1);
        if ~isempty(lagCI) && ~any(isnan(lagCI))
            hLag1 = xline(lagCI(1), '--', 'Color', lagCIColor, 'LineWidth', 1.4);
            xline(lagCI(2), '--', 'Color', lagCIColor, 'LineWidth', 1.4);
        end

        xlabel('Lag (Seconds)');
        ylabel('Correlation');
        title(sprintf('%s\nbouts = %d | dur = %.1f s', ...
            goodBehLabels{k}, goodNBouts(k), goodDurSec(k)));
        box off;
        set(gca, 'FontSize', 14, 'LineWidth', 1, 'TickDir', 'out');

        if k == 1
            legHandles = [hOrig];
            legLabels = {'Original Cross-Correlation'};
            if isgraphics(hCorr25), legHandles(end+1) = hCorr25; legLabels{end+1} = '2.5% Correlation Control'; end
            if isgraphics(hCorr97), legHandles(end+1) = hCorr97; legLabels{end+1} = '97.5% Correlation Control'; end
            legHandles(end+1) = hPeak; legLabels{end+1} = 'Actual Peak Lag';
            if isgraphics(hLag1), legHandles(end+1) = hLag1; legLabels{end+1} = 'Lag 95% CI (Permutation)'; end

            lgd = legend(legHandles, legLabels, 'Orientation', 'horizontal');
            lgd.Layout.Tile = 'south';
            lgd.FontSize = 12;
            lgd.Box = 'off';
        end
    end

    title(tile_lay, sprintf('%s nonchunked: lag vs correlation by classifier behavior', animalIDs{s}), ...
        'FontSize', 16, 'FontWeight', 'bold');

    % ---- figure 2: permutation histograms ----
    figure('Name', sprintf('%s nonchunk permutation histograms by behavior', animalIDs{s}), 'Color', 'w');
    tile_lay = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

    for k = 1:nPlot
        nexttile(tile_lay, k); hold on;

        permLags = goodPermLags{k};
        thisBeh = goodBehStruct{k};

        if isempty(permLags)
            title(sprintf('%s\n(no perms)', goodBehLabels{k}));
            axis off;
            continue
        end

        permLags = permLags(~isnan(permLags));
        histogram(permLags, 'FaceColor', [0.3 0.6 0.8], 'EdgeColor', 'none');

        prcLag = prctile(permLags, [2.5 97.5]);
        hCI1 = xline(prcLag(1), '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
        xline(prcLag(2), '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
        hPeak = xline(thisBeh.peakLagSec, 'r-', 'LineWidth', 1.5);

        xlabel('Peak Lag (s)');
        ylabel('Count');
        title(sprintf('%s\nbouts = %d | dur = %.1f s', ...
            goodBehLabels{k}, goodNBouts(k), goodDurSec(k)));
        box off;
        set(gca, 'FontSize', 14, 'LineWidth', 1, 'TickDir', 'out');

        if k == 1
            lgd = legend([hCI1 hPeak], {'Permutation 95% CI', 'Actual Peak Lag'}, ...
                'Orientation', 'horizontal');
            lgd.Layout.Tile = 'south';
            lgd.FontSize = 12;
            lgd.Box = 'off';
        end
    end

    title(tile_lay, sprintf('%s nonchunked: permutation lag distributions by classifier behavior', animalIDs{s}), ...
        'FontSize', 16, 'FontWeight', 'bold');
end

end
