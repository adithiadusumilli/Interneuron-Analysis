function plotSavedChunkedPopavgCrossCorrelationResults(saveFile)
% plots saved chunked population-averaged cross-correlation outputs
% without rerunning the analysis

% j run: plotSavedChunkedPopavgCrossCorrelationResults("C:\Users\mirilab\Documents\GlobusTransfer\concatCrossCorr_chunked_popavg_ALL_SESSIONS.mat")

arguments
    saveFile (1,1) string
end

S = load(saveFile, 'allSessions');
allSessions = S.allSessions;

sessions = allSessions.sessions;
nSess = numel(sessions);

% extract animal IDs
animalIDs = cell(1,nSess);
for i = 1:nSess
    a = regexp(sessions(i).baseDir,'D\d+','match','once');
    if isempty(a)
        a = sprintf('Session %d',i);
    end
    animalIDs{i} = a;
end

% colors
origColor = [0 0 0];
corrCIColor = [0 0.2 0.6];
peakLagColor = [0.95 0.45 0.35];
lagCIColor = [0 0 0];

peakLags = nan(nSess,1);
lagCIAll = nan(nSess,2);
permLagCell = cell(nSess,1);

%% ---- tiled cross-correlation figure ----
figure('Name', 'M1 Lag vs. Correlation', 'Color', 'w');
tile_lay = tiledlayout(1, nSess, 'TileSpacing', 'compact', 'Padding', 'compact');

for iDir = 1:nSess
    nexttile(tile_lay, iDir); hold on;

    sess = sessions(iDir);

    if isempty(sess.real_xc) || isempty(sess.lags)
        title(sprintf('%s (missing)', animalIDs{iDir}));
        axis off;
        continue;
    end

    lagsSec = sess.lags * sess.binSize;
    xc = sess.real_xc(:)';
    peakLag = sess.real_peakLagSec;
    peakLags(iDir) = peakLag;

    if ~isempty(sess.shift_xcZeroLag)
        corrCI = prctile(sess.shift_xcZeroLag, [2.5 97.5]);
    else
        corrCI = [NaN NaN];
    end

    if ~isempty(sess.perm_peakLagSec)
        goodPerm = ~isnan(sess.perm_peakLagSec);
        if any(goodPerm)
            lagCI = prctile(sess.perm_peakLagSec(goodPerm), [2.5 97.5]);
            permLagCell{iDir} = sess.perm_peakLagSec(goodPerm);
        else
            lagCI = [NaN NaN];
        end
    else
        lagCI = [NaN NaN];
    end
    lagCIAll(iDir,:) = lagCI;

    hOrig = plot(lagsSec, xc, 'Color', origColor, 'LineWidth', 2);

    if ~any(isnan(corrCI))
        hCorr25 = yline(corrCI(1), '--', 'Color', corrCIColor, 'LineWidth', 1.4);
        hCorr97 = yline(corrCI(2), '--', 'Color', corrCIColor, 'LineWidth', 1.4);
    end

    hPeak = xline(peakLag, '--', 'Color', peakLagColor, 'LineWidth', 1.8);

    if ~any(isnan(lagCI))
        hLag1 = xline(lagCI(1), '--', 'Color', lagCIColor, 'LineWidth', 1.4);
        hLag2 = xline(lagCI(2), '--', 'Color', lagCIColor, 'LineWidth', 1.4);
    end

    xlabel('Lag (Seconds)');
    ylabel('Correlation');
    title(sprintf('M1 Lag vs. Correlation — %s', animalIDs{iDir}));
    box off;
    set(gca, 'FontSize', 18, 'LineWidth', 1, 'TickDir', 'out');

    if iDir == 1
        lgd = legend([hOrig hCorr25 hCorr97 hPeak hLag1], ...
            {'Original Cross-Correlation', ...
             '2.5% Correlation Control', ...
             '97.5% Correlation Control', ...
             'Actual Peak Lag', ...
             'Lag 95% CI (Permutation)'}, ...
            'Orientation', 'horizontal');

        lgd.Layout.Tile = 'south';
        lgd.FontSize = 16;
        lgd.Box = 'off';
    end
end

%% ---- summary figure: peak lag + ci per animal ----
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
xlabel('Animal');
ylabel('Peak Lag (Seconds)');
xticks(xPos);
xticklabels(animalIDs);

title('M1 Peak Lags with 95% Permutation CI');
box off;
grid on;

%% ---- per-animal permutation histograms ----
figure('Name','M1 Chunked XC Peak Lag Permutation Distributions (Per Animal)', ...
       'Color','w');
tiledlayout(1, nSess, 'TileSpacing','compact','Padding','compact');

for s = 1:nSess
    nexttile; hold on;

    permLags = permLagCell{s};

    if isempty(permLags)
        title(sprintf('%s (no perms)', animalIDs{s}));
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

    title(animalIDs{s});
    box off;
end

%% ---- combined histogram across animals ----
allPermLags = cat(2, permLagCell{:});
validAll = ~isnan(allPermLags);
allPermLags = allPermLags(validAll);

figure('Name','M1 Chunked XC Peak Lag Permutations (All Animals Combined)', ...
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
        actualLines(end+1,1) = xline(peakLags(s), '-', ...
            'Color',co(s,:), 'LineWidth',1.5); %#ok<AGROW>
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
    legEntries{end+1} = sprintf('%s Lag = %.3g s', animalIDs{s}, peakLags(s));
end

legend([hHist h2_5 h97_5 actualLines(:)'], legEntries, 'Location', 'best');

end
