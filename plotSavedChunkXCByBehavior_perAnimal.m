function plotSavedChunkXCByBehavior_perAnimal(chunkFile, minTrials)
% plots saved CHUNKED popavg xc results
% split by classifier behv

% for each animal, makes 2 figures:
%   1. lag vs correlation, with 1 panel per behv that passed minTrials
%   2. permutation lag histograms, with 1 panel per behv that passed minTrials

% j run: plotSavedChunkXCByBehavior_perAnimal("C:\Users\mirilab\Documents\GlobusTransfer\combined_allAnimals_concatCrossCorrPerCanonicalBehavior_classifier.mat", 100)

arguments
    chunkFile (1,1) string
    minTrials (1,1) double = 100
end

behNames = {'unlabeled','climbdown','climbup','eating','grooming','jumpdown','jumping','rearing','still','walkflat','walkgrid'};

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

% force animal ids to be D026 / D020 / D024 / D043 from baseDirs
for i = 1:nSess
    a = regexp(C.baseDirs{i}, 'D\d+', 'match', 'once');
    if isempty(a)
        a = sprintf('Session %d', i);
    end
    animalIDs{i} = a;
end

fprintf('\n================ chunked per animal ================\n')

for s = 1:nSess
    goodBehLabels = {};
    goodXC = {};
    goodPeakLags = [];
    goodCorrCI = {};
    goodLagCI = {};
    goodNTrials = [];
    goodPermLags = {};

    for bIdx = 1:numel(C.behaviors)
        beh = C.behaviors(bIdx);

        nTrials = C.all_nTrials_real{s, bIdx};
        if isempty(nTrials) || isnan(nTrials) || nTrials < minTrials
            continue
        end

        xc = C.all_xc_real{s, bIdx};
        if isempty(xc)
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

        if beh >= 0 && beh <= 10
            goodBehLabels{end+1} = sprintf('%d: %s', beh, behNames{beh+1}); %#ok<AGROW>
        else
            goodBehLabels{end+1} = sprintf('%d', beh); %#ok<AGROW>
        end

        goodXC{end+1} = xc; %#ok<AGROW>
        goodPeakLags(end+1) = C.all_peakLagSec_real(s, bIdx); %#ok<AGROW>
        goodCorrCI{end+1} = thisCorrCI; %#ok<AGROW>
        goodLagCI{end+1} = thisLagCI; %#ok<AGROW>
        goodNTrials(end+1) = nTrials; %#ok<AGROW>
        goodPermLags{end+1} = permLags; %#ok<AGROW>
    end

    fprintf('%s chunked: %d behaviors passed minTrials=%d\n', ...
        animalIDs{s}, numel(goodBehLabels), minTrials)

    if isempty(goodBehLabels)
        continue
    end

    nPlot = numel(goodBehLabels);
    nCols = min(4, nPlot);
    nRows = ceil(nPlot / nCols);

    % ---- figure 1: lag vs correlation ----
    figure('Name', sprintf('%s chunked lag vs correlation by behavior', animalIDs{s}), 'Color', 'w');
    tile_lay = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

    for k = 1:nPlot
        nexttile(tile_lay, k); hold on;

        lagsSec = C.lags * C.binSize;
        xc = goodXC{k};
        peakLag = goodPeakLags(k);
        corrCI = goodCorrCI{k};
        lagCI = goodLagCI{k};

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
        title(sprintf('%s\nn = %d', goodBehLabels{k}, goodNTrials(k)));
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

    title(tile_lay, sprintf('%s chunked: lag vs correlation by classifier behavior', animalIDs{s}), ...
        'FontSize', 16, 'FontWeight', 'bold');

    % ---- figure 2: permutation histograms ----
    figure('Name', sprintf('%s chunked permutation histograms by behavior', animalIDs{s}), 'Color', 'w');
    tile_lay = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

    for k = 1:nPlot
        nexttile(tile_lay, k); hold on;

        permLags = goodPermLags{k};

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
        hPeak = xline(goodPeakLags(k), 'r-', 'LineWidth', 1.5);

        xlabel('Peak Lag (s)');
        ylabel('Count');
        title(sprintf('%s\nn = %d', goodBehLabels{k}, goodNTrials(k)));
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

    title(tile_lay, sprintf('%s chunked: permutation lag distributions by classifier behavior', animalIDs{s}), ...
        'FontSize', 16, 'FontWeight', 'bold');
end

end
