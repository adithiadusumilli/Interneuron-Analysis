function plotSavedChunkXC_SelectedBehaviors_AllAnimalsSummary(chunkFile, minTrials)
% plots saved chunked popavg xc results for selected classifier behaviors: climbup (2), eating (3), walkgrid (10)

% outputs:
%   1. per animal: 1x3 lag vs corr plots
%   2. per animal: 1x3 permutation lag histograms
%   3. one summary 1x3 figure:
%        each panel = one behavior
%        each panel contains all 4 animals:
%           - black dot = actual peak lag
%           - slightly offset gray vertical line = null 95% range

% j run: plotSavedChunkXC_SelectedBehaviors_AllAnimalsSummary("C:\Users\mirilab\Documents\GlobusTransfer\combined_allAnimals_concatCrossCorrPerCanonicalBehavior_classifier.mat", 25)

arguments
    chunkFile (1,1) string
    minTrials (1,1) double = 25
end

behNums = [2 3 10];
behNames = {'climbup','eating','walkgrid'};

origColor = [0 0 0];
corrCIColor = [0 0.2 0.6];
peakLagColor = [0.95 0.45 0.35];
lagCIColor = [0 0 0];

C = load(chunkFile, ...
    'animalNames', 'baseDirs', 'behaviors', 'lags', 'binSize', ...
    'all_xc_real', 'all_peakLagSec_real', 'all_nTrials_real', ...
    'all_peakLagSec_perm', 'all_xcZeroLag_shift');

nSess = numel(C.animalNames);
animalIDs = cell(1, nSess);

for i = 1:nSess
    a = regexp(C.baseDirs{i}, 'D\d+', 'match', 'once');
    if isempty(a)
        a = sprintf('Session %d', i);
    end
    animalIDs{i} = a;
end

actualLagMat = nan(nSess, numel(behNums));
lagCIMatLo = nan(nSess, numel(behNums));
lagCIMatHi = nan(nSess, numel(behNums));
trialMat = nan(nSess, numel(behNums));

fprintf('\n================ chunked selected behavior plots ================\n')

for s = 1:nSess
    % ---- lag vs correlation ----
    figure('Name', sprintf('%s chunked lag vs corr selected behaviors', animalIDs{s}), 'Color', 'w');
    tile_lay1 = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    firstLegendDone = false;

    for k = 1:numel(behNums)
        b = behNums(k);
        bIdx = find(C.behaviors == b, 1);

        nexttile(tile_lay1, k); hold on;

        if isempty(bIdx)
            title(sprintf('%d: %s\nmissing', b, behNames{k}));
            axis off;
            continue
        end

        nTrials = C.all_nTrials_real(s, bIdx);
        xc = C.all_xc_real{s, bIdx};

        if isempty(nTrials) || isnan(nTrials) || nTrials < minTrials || isempty(xc)
            title(sprintf('%d: %s\nn = %g (< %d)', b, behNames{k}, nTrials, minTrials));
            axis off;
            continue
        end

        permLags = C.all_peakLagSec_perm{s, bIdx};
        shiftZero = C.all_xcZeroLag_shift{s, bIdx};

        thisCorrCI = [NaN NaN];
        if ~isempty(shiftZero)
            shiftZero = shiftZero(~isnan(shiftZero));
            if ~isempty(shiftZero)
                thisCorrCI = prctile(shiftZero, [2.5 97.5]);
            end
        end

        thisLagCI = [NaN NaN];
        if ~isempty(permLags)
            permLags = permLags(~isnan(permLags));
            if ~isempty(permLags)
                thisLagCI = prctile(permLags, [2.5 97.5]);
            end
        end

        lagsSec = C.lags * C.binSize;
        peakLag = C.all_peakLagSec_real(s, bIdx);

        actualLagMat(s,k) = peakLag;
        if ~any(isnan(thisLagCI))
            lagCIMatLo(s,k) = thisLagCI(1);
            lagCIMatHi(s,k) = thisLagCI(2);
        end
        trialMat(s,k) = nTrials;

        hOrig = plot(lagsSec, xc, 'Color', origColor, 'LineWidth', 2);

        hCorr25 = gobjects(1);
        hCorr97 = gobjects(1);
        if ~any(isnan(thisCorrCI))
            hCorr25 = yline(thisCorrCI(1), '--', 'Color', corrCIColor, 'LineWidth', 1.4);
            hCorr97 = yline(thisCorrCI(2), '--', 'Color', corrCIColor, 'LineWidth', 1.4);
        end

        hPeak = xline(peakLag, '--', 'Color', peakLagColor, 'LineWidth', 1.8);

        hLag1 = gobjects(1);
        if ~any(isnan(thisLagCI))
            hLag1 = xline(thisLagCI(1), '--', 'Color', lagCIColor, 'LineWidth', 1.4);
            xline(thisLagCI(2), '--', 'Color', lagCIColor, 'LineWidth', 1.4);
        end

        xlabel('Lag (Seconds)');
        ylabel('Correlation');
        title(sprintf('%d: %s\nn = %d', b, behNames{k}, nTrials));
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

    title(tile_lay1, sprintf('%s chunked: selected behaviors lag vs correlation', animalIDs{s}), ...
        'FontSize', 16, 'FontWeight', 'bold');

    % ---- permutation lag histograms ----
    figure('Name', sprintf('%s chunked permutation selected behaviors', animalIDs{s}), 'Color', 'w');
    tile_lay2 = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    firstLegendDone = false;

    for k = 1:numel(behNums)
        b = behNums(k);
        bIdx = find(C.behaviors == b, 1);

        nexttile(tile_lay2, k); hold on;

        if isempty(bIdx)
            title(sprintf('%d: %s\nmissing', b, behNames{k}));
            axis off;
            continue
        end

        nTrials = C.all_nTrials_real(s, bIdx);
        permLags = C.all_peakLagSec_perm{s, bIdx};

        if isempty(nTrials) || isnan(nTrials) || nTrials < minTrials || isempty(permLags)
            title(sprintf('%d: %s\nn = %g (< %d)', b, behNames{k}, nTrials, minTrials));
            axis off;
            continue
        end

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
        hPeak = xline(C.all_peakLagSec_real(s, bIdx), 'r-', 'LineWidth', 1.5);

        xlabel('Peak Lag (s)');
        ylabel('Count');
        title(sprintf('%d: %s\nn = %d', b, behNames{k}, nTrials));
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

    title(tile_lay2, sprintf('%s chunked: selected behavior permutation lag distributions', animalIDs{s}), ...
        'FontSize', 16, 'FontWeight', 'bold');
end

% ---- summary figure ----
figure('Name', 'Chunked selected behavior summary', 'Color', 'w');
tile_lay3 = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

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

title(tile_lay3, 'Chunked summary: actual peak lag with offset null 95% range by selected behavior', ...
    'FontSize', 16, 'FontWeight', 'bold');

end
