function runCrossCorrelation_plotOnly_getMouseDataNames(saveFile)
% plot-only function for runCrossCorrelation_saveOnly_getMouseDataNames
% loads saved cross-correlation outputs and regenerates summary + permutation plots

if nargin < 1 || isempty(saveFile)
    saveFile = fullfile('X:\David\AnalysesData\InterneuronAnalyses\Lab Meeting Pres', ...
        '6 animals run cross correlation, no chunking, pop-wise', ...
        'runCrossCorrelation_savedOutputs_all6Animals.mat');
end

S = load(saveFile);

xcorrResults = S.xcorrResults;
peakLags = S.peakLags;
peakCorrs = S.peakCorrs;
lagCIAll = S.lagCIAll;
permLagCell = S.permLagCell;

nSess = numel(peakLags);

if isfield(S, 'animalLabels')
    animalLabels = S.animalLabels;
else
    animalLabels = cell(1,nSess);
    for s = 1:nSess
        animalLabels{s} = sprintf('Session %d', s);
    end
end

if iscell(animalLabels)
    animalLabels = string(animalLabels);
end

%% ============================
%  plot cross-correlation curves
% ============================

figure('Name','Cortex population cross-correlation curves', 'Color','w');
tiledlayout(1, nSess, 'TileSpacing','compact','Padding','compact');

for s = 1:nSess
    nexttile; hold on;

    sess = xcorrResults.sessions(s);

    if isempty(sess.xc) || isempty(sess.lagsSec)
        title(sprintf('%s missing', animalLabels(s)));
        axis off;
        continue;
    end

    plot(sess.lagsSec, sess.xc, 'k', 'LineWidth', 2);
    xline(sess.peakLag, 'r-', 'LineWidth', 1.5);
    yline(sess.corrCI(1), '--', 'Color',[0.2 0.2 0.2], 'LineWidth',1);
    yline(sess.corrCI(2), '--', 'Color',[0.2 0.2 0.2], 'LineWidth',1);

    xlabel('Lag (s)');
    ylabel('Correlation');
    title(sprintf('%s', animalLabels(s)));
    box off;
end

%% ============================
%  summary peak lag + CI
% ============================

figure('Name','Cortex peak lag summary with permutation CI', 'Color', 'w'); 
hold on;

xPos = 1:nSess;

for i = 1:nSess
    if ~any(isnan(lagCIAll(i,:)))
        line([xPos(i) xPos(i)], lagCIAll(i,:), ...
            'Color', [0.6 0.6 0.6], 'LineWidth', 2);
    end
end

scatter(xPos, peakLags, 70, 'k', 'filled');
yline(0, 'k:');

xlim([0.5 nSess + 0.5]);
xlabel('Animal');
ylabel('Peak lag (seconds)');
xticks(xPos);
xticklabels(cellstr(animalLabels));

title('Cortex peak lags with 95% permutation CI');
box off;
grid on;

fprintf('\n========== summary cortex ==========\n');
fprintf('peak lags (s):\n');  disp(peakLags);
fprintf('peak corrs:\n');     disp(peakCorrs);
fprintf('lag CIs (s):\n');    disp(lagCIAll);

%% ============================
%  common histogram limits
% ============================

allPermForLimits = cat(2, permLagCell{:});
allPermForLimits = allPermForLimits(~isnan(allPermForLimits));

allActualForLimits = peakLags(~isnan(peakLags));

if isempty(allPermForLimits) && isempty(allActualForLimits)
    commonXLim = [-0.5 0.5];
else
    allVals = [allPermForLimits(:); allActualForLimits(:)];
    xMin = min(allVals);
    xMax = max(allVals);

    if xMin == xMax
        pad = 0.001;
    else
        pad = 0.05 * (xMax - xMin);
    end

    commonXLim = [xMin - pad, xMax + pad];
end

nBins = 24;
commonEdges = linspace(commonXLim(1), commonXLim(2), nBins + 1);

%% ============================
%  per-animal permutation histograms
% ============================

figure('Name','Cortex XC peak lag permutation distributions per animal', ...
       'Color','w');

tiledlayout(1, nSess, 'TileSpacing','compact','Padding','compact');

for s = 1:nSess
    nexttile; hold on;

    permLags = permLagCell{s};
    permLags = permLags(~isnan(permLags));

    if isempty(permLags)
        title(sprintf('%s no perms', animalLabels(s)));
        axis off;
        continue;
    end

    histogram(permLags, ...
        'BinEdges', commonEdges, ...
        'FaceColor',[0.3 0.6 0.8], ...
        'EdgeColor','none');

    prcLag = prctile(permLags, [2.5 97.5]);

    xline(prcLag(1), '--', 'Color',[0.2 0.2 0.2], 'LineWidth',1.5);
    xline(prcLag(2), '--', 'Color',[0.2 0.2 0.2], 'LineWidth',1.5);
    xline(peakLags(s), 'r-', 'LineWidth',1.5);

    xlim(commonXLim);
    xlabel('Peak lag (s)');
    ylabel('Count');
    title(sprintf('%s', animalLabels(s)));
    box off;
end

%% ============================
%  combined histogram
% ============================

allPermLags = cat(2, permLagCell{:});
allPermLags = allPermLags(~isnan(allPermLags));

figure('Name','Cortex XC peak lag permutations all animals combined', ...
       'Color','w'); 
hold on;

hHist = histogram(allPermLags, ...
    'BinEdges', commonEdges, ...
    'FaceColor',[0.3 0.6 0.8], ...
    'EdgeColor','none');

xlabel('Peak lag (s)');
ylabel('Count');

prcAll = prctile(allPermLags, [2.5 97.5]);
h2_5 = xline(prcAll(1), '--', 'Color',[0.2 0.2 0.2], 'LineWidth',1.5);
h97_5 = xline(prcAll(2), '--', 'Color',[0.2 0.2 0.2], 'LineWidth',1.5);

co = lines(nSess);
actualLines = gobjects(0,1);

for s = 1:nSess
    if ~isnan(peakLags(s))
        actualLines(end+1,1) = xline(peakLags(s), '-', ...
            'Color', co(s,:), 'LineWidth',1.5); %#ok<AGROW>
    end
end

xlim(commonXLim);
title('Combined permutation distribution across all animals');
box off;

legHandles = [hHist h2_5 h97_5 actualLines(:)'];

legEntries = { ...
    'Permuted lags', ...
    sprintf('2.5%% all = %.3g s', prcAll(1)), ...
    sprintf('97.5%% all = %.3g s', prcAll(2))};

for s = 1:nSess
    legEntries{end+1} = sprintf('%s actual = %.3g s', animalLabels(s), peakLags(s)); %#ok<AGROW>
end

legend(legHandles, legEntries, 'Location','best', 'Box','off');

end
