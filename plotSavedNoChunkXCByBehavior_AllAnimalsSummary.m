function plotSavedNoChunkXCByBehavior_AllAnimalsSummary(noChunkFile)
% plots saved nonchunked popavg xc results by behavior

% outputs:
%   1. per animal: 2x5 lag vs corr plots (behaviors 1:10)
%   2. per animal: 2x5 permutation lag histograms (behaviors 1:10)
%   3. one summary 2x5 figure:
%        each panel = one behavior
%        each panel contains all 4 animals:
%           - black dot = actual peak lag
%           - slightly offset gray vertical line = null 95% range

% this version plots behaviors regardless of bout count, as long as data exist

% j run: plotSavedNoChunkXCByBehavior_AllAnimalsSummary("X:\David\AnalysesData\nonchunked_xcorr_by_classifier_cortex_allSessions_saved.mat")

arguments
    noChunkFile (1,1) string
end

behNums = 1:10;
behNames = {'climbdown','climbup','eating','grooming','jumpdown', ...
            'jumping','rearing','still','walkflat','walkgrid'};

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

binSize = 0.001;

% summary storage: session x behavior
actualLagMat = nan(nSess, numel(behNums));
lagCIMatLo = nan(nSess, numel(behNums));
lagCIMatHi = nan(nSess, numel(behNums));
boutMat = nan(nSess, numel(behNums));
durMat = nan(nSess, numel(behNums));

fprintf('\n================ nonchunked behavior plots ================\n')

for s = 1:nSess
    figure('Name', sprintf('%s nonchunk lag vs corr by behavior', animalIDs{s}), 'Color', 'w');
    tile_lay1 = tiledlayout(2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

    firstLegendDone = false;

    for k = 1:numel(behNums)
        b = behNums(k);
        thisBeh = R.sessions(s).beh(b);

        nexttile(tile_lay1, k); hold on;

        if isempty(thisBeh) || ~isfield(thisBeh, 'timeIdx') || isempty(thisBeh.timeIdx) || ...
                isempty(thisBeh.xc) || isempty(thisBeh.lagsSec)
            title(sprintf('%d: %s\nmissing', b, behNames{k}));
            axis off;
            continue
        end

        timeIdx = sort(thisBeh.timeIdx(:));
        if isempty(timeIdx)
            nBouts = 0;
        else
            nBouts = 1 + sum(diff(timeIdx) > 1);
        end
        durSec = numel(timeIdx) * binSize;

        lagsSec = thisBeh.lagsSec(:)';
        xc = thisBeh.xc(:)';
        peakLag = thisBeh.peakLagSec;
        corrCI = thisBeh.ctrlCorrCI;
        lagCI = thisBeh.lagCI;

        actualLagMat(s,k) = peakLag;
        if ~isempty(lagCI) && numel(lagCI) == 2 && ~any(isnan(lagCI))
            lagCIMatLo(s,k) = lagCI(1);
            lagCIMatHi(s,k) = lagCI(2);
        end
        boutMat(s,k) = nBouts;
        durMat(s,k) = durSec;

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
        title(sprintf('%d: %s\nbouts = %d | dur = %.1f s', ...
            b, behNames{k}, nBouts, durSec));
        box off;
        set(gca, 'FontSize', 12, 'LineWidth', 1, 'TickDir', 'out');

        if ~firstLegendDone
            legHandles = [hOrig];
            legLabels = {'Original Cross-Correlation'};
            if isgraphics(hCorr25), legHandles(end+1) = hCorr25; legLabels{end+1} = '2.5% Correlation Control'; end
            if isgraphics(hCorr97), legHandles(end+1) = hCorr97; legLabels{end+1} = '97.5% Correlation Control'; end
            legHandles(end+1) = hPeak; legLabels{end+1} = 'Actual Peak Lag';
            if isgraphics(hLag1), legHandles(end+1) = hLag1; legLabels{end+1} = 'Lag 95% CI (Permutation)'; end

            lgd = legend(legHandles, legLabels, 'Orientation', 'horizontal');
            lgd.Layout.Tile = 'south';
            lgd.FontSize = 11;
            lgd.Box = 'off';
            firstLegendDone = true;
        end
    end

    title(tile_lay1, sprintf('%s nonchunked: lag vs correlation by behavior', animalIDs{s}), ...
        'FontSize', 16, 'FontWeight', 'bold');

    % ---- permutation histograms ----
    figure('Name', sprintf('%s nonchunk permutation histograms by behavior', animalIDs{s}), 'Color', 'w');
    tile_lay2 = tiledlayout(2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

    firstLegendDone = false;

    for k = 1:numel(behNums)
        b = behNums(k);
        thisBeh = R.sessions(s).beh(b);

        nexttile(tile_lay2, k); hold on;

        if isempty(thisBeh) || ~isfield(thisBeh, 'timeIdx') || isempty(thisBeh.timeIdx) || ...
                isempty(thisBeh.permPeakLags)
            title(sprintf('%d: %s\nmissing', b, behNames{k}));
            axis off;
            continue
        end

        timeIdx = sort(thisBeh.timeIdx(:));
        if isempty(timeIdx)
            nBouts = 0;
        else
            nBouts = 1 + sum(diff(timeIdx) > 1);
        end
        durSec = numel(timeIdx) * binSize;

        permLags = thisBeh.permPeakLags;
        permLags = permLags(~isnan(permLags));

        if isempty(permLags)
            title(sprintf('%d: %s\n(no perms)', b, behNames{k}));
            axis off;
            continue
        end

        hHist = histogram(permLags, 'FaceColor', [0.3 0.6 0.8], 'EdgeColor', 'none');

        prcLag = prctile(permLags, [2.5 97.5]);
        hCI1 = xline(prcLag(1), '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
        xline(prcLag(2), '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
        hPeak = xline(thisBeh.peakLagSec, 'r-', 'LineWidth', 1.5);

        xlabel('Peak Lag (s)');
        ylabel('Count');
        title(sprintf('%d: %s\nbouts = %d | dur = %.1f s', ...
            b, behNames{k}, nBouts, durSec));
        box off;
        set(gca, 'FontSize', 12, 'LineWidth', 1, 'TickDir', 'out');

        if ~firstLegendDone
            lgd = legend([hHist hCI1 hPeak], ...
                {'Permutation Lag Distribution', 'Permutation 95% CI', 'Actual Peak Lag'}, ...
                'Orientation', 'horizontal');
            lgd.Layout.Tile = 'south';
            lgd.FontSize = 11;
            lgd.Box = 'off';
            firstLegendDone = true;
        end
    end

    title(tile_lay2, sprintf('%s nonchunked: permutation lag distributions by behavior', animalIDs{s}), ...
        'FontSize', 16, 'FontWeight', 'bold');
end

% ---- summary figure: one panel per behavior, all animals in each panel ----
figure('Name', 'Nonchunked summary by behavior', 'Color', 'w');
tile_lay3 = tiledlayout(2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

xBase = 1:nSess;
xOffsets = [-0.18 -0.06 0.06 0.18];

for k = 1:numel(behNums)
    nexttile(tile_lay3, k); hold on;

    for s = 1:nSess
        if ~isnan(lagCIMatLo(s,k)) && ~isnan(lagCIMatHi(s,k))
            xci = xBase(s) + xOffsets(s);
            line([xci xci], [lagCIMatLo(s,k) lagCIMatHi(s,k)], ...
                'Color', [0.6 0.6 0.6], 'LineWidth', 2);
        end

        if ~isnan(actualLagMat(s,k))
            plot(xBase(s), actualLagMat(s,k), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
        end
    end

    yline(0, 'k:');
    xlim([0.5 nSess + 0.5]);
    xticks(1:nSess);
    xticklabels(animalIDs);
    ylabel('Peak Lag (s)');
    title(sprintf('%d: %s', behNums(k), behNames{k}));
    box off
    grid on
end

hActual = plot(nan, nan, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
hNull = line([nan nan], [nan nan], 'Color', [0.6 0.6 0.6], 'LineWidth', 2);

lgd = legend([hActual hNull], ...
    {'Actual peak lag', 'Null 95% range (not CI on dot)'}, ...
    'Orientation', 'horizontal');
lgd.Layout.Tile = 'south';
lgd.FontSize = 11;
lgd.Box = 'off';

title(tile_lay3, 'Nonchunked summary: actual peak lag with offset null 95% range by behavior', ...
    'FontSize', 16, 'FontWeight', 'bold');

end
