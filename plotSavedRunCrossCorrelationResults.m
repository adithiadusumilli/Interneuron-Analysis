function plotSavedRunCrossCorrelationResults(saveFile)
% plotSavedCrossCorrelationResults

% loads saved run cross-correlation outputs (no chunking, pop-avgd) and redraws figures without rerunning analysis

% j run:
% plotSavedCrossCorrelationResults("X:\David\AnalysesData\InterneuronAnalyses\Lab Meeting Pres\4 aninals run cross correlation, no chunking, pop-wise\runCrossCorrelation_savedOutputs_all4Animals.mat")

arguments
    saveFile (1,1) string
end

S = load(saveFile, 'xcorrResults', 'peakLags', 'peakCorrs', 'lagCIAll', 'permLagCell', 'animalLabels');

xcorrResults = S.xcorrResults;
peakLags = S.peakLags;
peakCorrs = S.peakCorrs; %#ok<NASGU>
lagCIAll = S.lagCIAll;
permLagCell = S.permLagCell;
animalLabels = S.animalLabels;

nSess = numel(xcorrResults.sessions);

% colors
origColor = [0 0 0];
corrCIColor = [0 0.2 0.6];
peakLagColor = [0.95 0.45 0.35];
lagCIColor = [0 0 0];

%% ---- tiled cross-correlation figure ----
figure('Name', 'Cross-Correlation Summary (Cortex only)', 'Color', 'w');
tile_lay = tiledlayout(1, nSess, 'TileSpacing', 'compact', 'Padding', 'compact');

for iDir = 1:nSess
    nexttile(tile_lay, iDir); hold on;

    sess = xcorrResults.sessions(iDir);

    if isempty(sess.lagsSec) || isempty(sess.xc)
        title(sprintf('Session %d (missing)', iDir));
        axis off;
        continue;
    end

    hOrig = plot(sess.lagsSec, sess.xc, 'Color', origColor, 'LineWidth', 2);
    hCorr25 = yline(sess.corrCI(1), '--', 'Color', corrCIColor, 'LineWidth', 1.4);
    hCorr97 = yline(sess.corrCI(2), '--', 'Color', corrCIColor, 'LineWidth', 1.4);
    hPeak = xline(sess.peakLag, '--', 'Color', peakLagColor, 'LineWidth', 1.8);

    hLagCI = gobjects(1,2);
    if ~any(isnan(sess.lagCI))
        hLagCI(1) = xline(sess.lagCI(1), '--', 'Color', lagCIColor, 'LineWidth', 1.4);
        hLagCI(2) = xline(sess.lagCI(2), '--', 'Color', lagCIColor, 'LineWidth', 1.4);
    end

    xlabel('Lag (Seconds)');
    ylabel('Correlation');
    title(sprintf('M1 Lag vs. Correlation — %s', sess.animalLabel));
    box off;
    set(gca, 'FontSize', 18, 'LineWidth', 1, 'TickDir', 'out');

    if iDir == 1
        lgd = legend([hOrig hCorr25 hCorr97 hPeak hLagCI(1)], ...
            {'Original Cross-Correlation', ...
             '2.5% Correlation Control', ...
             '97.5% Correlation Control', ...
             'Actual Peak Lag', ...
             'Lag 95% CI (Permutation)'}, ...
            'Orientation', 'horizontal');
        lgd.Layout.Tile = 'south';
        lgd.FontSize = 12;
        lgd.Box = 'off';
    end
end

%% ---- summary figure: peak lag + ci per session ----
figure('Name','Cortex peak lag summary with permutation CI', 'Color', 'w'); hold on;

xPos = 1:nSess;

for i = 1:nSess
    if ~any(isnan(lagCIAll(i,:)))
        line([xPos(i) xPos(i)], lagCIAll(i,:), 'Color', [0.6 0.6 0.6], 'LineWidth', 2);
    end
end

scatter(xPos, peakLags, 70, 'k', 'filled');
yline(0, 'k:');

xlim([0.5 nSess + 0.5]);
xlabel('Session');
ylabel('Peak Lag (Seconds)');
xticks(xPos);
xticklabels(arrayfun(@(k) sprintf('Session %d', k), 1:nSess, 'UniformOutput', false));

title('Cortex Peak Lags with 95% Permutation CI');
box off;
grid on;

%% ---- per-animal histograms ----
figure('Name','Cortex chunked XC peak lag permutation distributions (per animal)', ...
       'Color','w');
tiledlayout(1, nSess, 'TileSpacing','compact','Padding','compact');

for s = 1:nSess
    nexttile; hold on;

    permLags = permLagCell{s};
    if isempty(permLags)
        title(sprintf('%s (no perms)', animalLabels{s}));
        axis off;
        continue;
    end

    histogram(permLags, 'FaceColor',[0.3 0.6 0.8], 'EdgeColor','none');
    xlabel('Peak Lag (s)');
    ylabel('Count');

    prcLag = prctile(permLags, [2.5 97.5]);
    xline(prcLag(1), '--', 'Color',[0.2 0.2 0.2], 'LineWidth',1.5);
    xline(prcLag(2), '--', 'Color',[0.2 0.2 0.2], 'LineWidth',1.5);

    xline(peakLags(s), 'r-', 'LineWidth',1.5);

    title(sprintf('Session %d', s));
    box off;
end

%% ---- combined histogram across animals ----
allPermLags = cat(2, permLagCell{:});
validAll = ~isnan(allPermLags);
allPermLags = allPermLags(validAll);

figure('Name','Cortex chunked XC peak lag permutations (all animals combined)', ...
       'Color','w'); hold on;

hHist = histogram(allPermLags, 'FaceColor',[0.3 0.6 0.8], 'EdgeColor','none');
xlabel('Peak Lag (s)');
ylabel('Count');

prcAll = prctile(allPermLags, [2.5 97.5]);
h2_5 = xline(prcAll(1), '--', 'Color',[0.2 0.2 0.2], 'LineWidth',1.5);
h97_5 = xline(prcAll(2), '--', 'Color',[0.2 0.2 0.2], 'LineWidth',1.5);

co = lines(nSess);
actualLines = gobjects(0,1);
for s = 1:nSess
    if ~isnan(peakLags(s))
        actualLines(end+1,1) = xline(peakLags(s), '-', 'Color',co(s,:), 'LineWidth',1.5); %#ok<AGROW>
    end
end

title('Combined Permutation Distribution (All Animals)');
box off;

legEntries = {
    'Permuted Lags', ...
    sprintf('2.5%% (All) = %.3g s', prcAll(1)), ...
    sprintf('97.5%% (All) = %.3g s', prcAll(2))
};

for s = 1:nSess
    legEntries{end+1} = sprintf('Actual Session %d Lag = %.3g s', s, peakLags(s));
end

legend([hHist h2_5 h97_5 actualLines(:)'], legEntries, 'Location', 'best');

end
