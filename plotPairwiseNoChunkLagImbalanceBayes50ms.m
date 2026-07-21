function plotPairwiseNoChunkLagImbalanceBayes50ms()
% Plot saved pairwise no-chunk lag-imbalance null distributions.

% Loads:
%   C:\Users\mirilab\Documents\GlobusTransfer\pairwiseNoChunkLagImbalanceBayes50ms_StoreyCorrThresh.mat

% For each animal, plots:
%   H0 = randomized labels with no lag shift
%   H+50 = randomized labels with +50 ms applied to oriented lags
%   H-50 = randomized labels with -50 ms applied to oriented lags

% The actual observed real int-pyr lag imbalance is shown as a red line.

% Also saves figures into:
%   C:\Users\mirilab\Documents\GlobusTransfer\
%   pairwiseLagImbalancePlots

%% ========================================================================
%% settings
%% ========================================================================

baseFolder = 'C:\Users\mirilab\Documents\GlobusTransfer';

inputFile = fullfile(baseFolder, 'pairwiseNoChunkLagImbalanceBayes50ms_StoreyCorrThresh.mat');

outputFolder = fullfile(baseFolder, 'pairwiseLagImbalancePlots');

if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

% Histogram settings
histEdges = -1:0.05:1;

% Set true to save individual and summary figures.
saveFigures = true;

% Set true to close each individual figure after saving.
closeIndividualFigures = false;

%% ========================================================================
%% load results
%% ========================================================================

if ~isfile(inputFile)
    error('Saved results file was not found:\n%s', inputFile);
end

loadedData = load(inputFile);

if ~isfield(loadedData, 'results')
    error('The MAT file does not contain a variable named results.');
end

results = loadedData.results;

if ~isfield(results, 'sessions') || isempty(results.sessions)
    error('results.sessions is missing or empty.');
end

fprintf('\nLoaded results from:\n%s\n', inputFile);

%% ========================================================================
%% identify valid sessions
%% ========================================================================

nSessionSlots = numel(results.sessions);
validSessionInds = [];

for sessInd = 1:nSessionSlots

    if isempty(results.sessions{sessInd})
        continue;
    end

    R = results.sessions{sessInd};

    requiredFields = { ...
        'animalID', ...
        'actualLagImbalance', ...
        'nullLagImbalanceH0', ...
        'nullLagImbalanceH50', ...
        'nullLagImbalanceHneg50'};

    hasAllFields = true;

    for fieldInd = 1:numel(requiredFields)

        if ~isfield(R, requiredFields{fieldInd})
            hasAllFields = false;
            break;
        end
    end

    if hasAllFields
        validSessionInds(end+1) = sessInd; %#ok<AGROW>
    end
end

if isempty(validSessionInds)
    error('No valid session results were found.');
end

nValidSessions = numel(validSessionInds);

fprintf('Found %d valid animal sessions.\n', nValidSessions);

%% ========================================================================
%% initialize summary values
%% ========================================================================

animalIDs = cell(nValidSessions,1);

actualImbalance = nan(nValidSessions,1);

medianH0 = nan(nValidSessions,1);
medianH50 = nan(nValidSessions,1);
medianHneg50 = nan(nValidSessions,1);

ciH0 = nan(nValidSessions,2);
ciH50 = nan(nValidSessions,2);
ciHneg50 = nan(nValidSessions,2);

pH0 = nan(nValidSessions,1);
pH50 = nan(nValidSessions,1);
pHneg50 = nan(nValidSessions,1);

ratioH0H50 = nan(nValidSessions,1);
ratioH0Hneg50 = nan(nValidSessions,1);

exactMatchRateH0 = nan(nValidSessions,1);
exactMatchRateH50 = nan(nValidSessions,1);
exactMatchRateHneg50 = nan(nValidSessions,1);

selectedCountMedianH0 = nan(nValidSessions,1);
selectedCountMedianH50 = nan(nValidSessions,1);
selectedCountMedianHneg50 = nan(nValidSessions,1);

%% ========================================================================
%% individual animal plots
%% ========================================================================

for validInd = 1:nValidSessions

    sessInd = validSessionInds(validInd);
    R = results.sessions{sessInd};

    animalID = char(string(R.animalID));
    animalIDs{validInd} = animalID;

    observed = R.actualLagImbalance;

    H0 = R.nullLagImbalanceH0(:);
    H50 = R.nullLagImbalanceH50(:);
    Hneg50 = R.nullLagImbalanceHneg50(:);

    H0 = H0(isfinite(H0));
    H50 = H50(isfinite(H50));
    Hneg50 = Hneg50(isfinite(Hneg50));

    actualImbalance(validInd) = observed;

    medianH0(validInd) = median(H0, 'omitnan');
    medianH50(validInd) = median(H50, 'omitnan');
    medianHneg50(validInd) = median(Hneg50, 'omitnan');

    if ~isempty(H0)
        ciH0(validInd,:) = prctile(H0, [2.5 97.5]);
    end

    if ~isempty(H50)
        ciH50(validInd,:) = prctile(H50, [2.5 97.5]);
    end

    if ~isempty(Hneg50)
        ciHneg50(validInd,:) = prctile(Hneg50, [2.5 97.5]);
    end

    if isfield(R, 'pValH0')
        pH0(validInd) = R.pValH0;
    end

    if isfield(R, 'pValH50')
        pH50(validInd) = R.pValH50;
    end

    if isfield(R, 'pValHneg50')
        pHneg50(validInd) = R.pValHneg50;
    end

    if isfield(R, 'evidenceRatio_H0_over_H50')
        ratioH0H50(validInd) = R.evidenceRatio_H0_over_H50;
    end

    if isfield(R, 'evidenceRatio_H0_over_Hneg50')
        ratioH0Hneg50(validInd) = R.evidenceRatio_H0_over_Hneg50;
    end

    %% ---------- exact-match information ----------

    if isfield(R, 'matchedTargetH0')
        exactMatchRateH0(validInd) = mean(R.matchedTargetH0, 'omitnan');
    end

    if isfield(R, 'matchedTargetH50')
        exactMatchRateH50(validInd) = mean(R.matchedTargetH50, 'omitnan');
    end

    if isfield(R, 'matchedTargetHneg50')
        exactMatchRateHneg50(validInd) = ...
            mean(R.matchedTargetHneg50, 'omitnan');
    end

    %% ---------- selected pair counts ----------

    if isfield(R, 'selectedLagCountsH0')
        selectedCountMedianH0(validInd) = ...
            median(R.selectedLagCountsH0, 'omitnan');
    end

    if isfield(R, 'selectedLagCountsH50')
        selectedCountMedianH50(validInd) = ...
            median(R.selectedLagCountsH50, 'omitnan');
    end

    if isfield(R, 'selectedLagCountsHneg50')
        selectedCountMedianHneg50(validInd) = ...
            median(R.selectedLagCountsHneg50, 'omitnan');
    end

    %% ====================================================================
    %% create individual figure
    %% ====================================================================

    fig = figure( ...
        'Name', sprintf('%s lag imbalance', animalID), ...
        'Color', 'w', ...
        'Position', [100 100 1500 460]);

    tiledlayout(1,3, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    %% ---------- H0 ----------

    nexttile;

    histogram(H0, histEdges, ...
        'Normalization', 'probability');

    hold on;

    xline(observed, 'r-', ...
        'LineWidth', 2.5, ...
        'Label', 'Actual', ...
        'LabelVerticalAlignment', 'bottom');

    xline(medianH0(validInd), 'k--', ...
        'LineWidth', 1.5, ...
        'Label', 'Null median');

    xlim([-1 1]);
    xlabel('Lag imbalance');
    ylabel('Proportion of permutations');

    title(sprintf( ...
        'H0: no shift\np = %.4f | median = %.3f', ...
        pH0(validInd), ...
        medianH0(validInd)));

    grid on;
    box off;

    %% ---------- H+50 ----------

    nexttile;

    histogram(H50, histEdges, ...
        'Normalization', 'probability');

    hold on;

    xline(observed, 'r-', ...
        'LineWidth', 2.5, ...
        'Label', 'Actual', ...
        'LabelVerticalAlignment', 'bottom');

    xline(medianH50(validInd), 'k--', ...
        'LineWidth', 1.5, ...
        'Label', 'Null median');

    xlim([-1 1]);
    xlabel('Lag imbalance');
    ylabel('Proportion of permutations');

    title(sprintf( ...
        'H+50 ms\np = %.4f | H0/H+50 = %.3f', ...
        pH50(validInd), ...
        ratioH0H50(validInd)));

    grid on;
    box off;

    %% ---------- H-50 ----------

    nexttile;

    histogram(Hneg50, histEdges, ...
        'Normalization', 'probability');

    hold on;

    xline(observed, 'r-', ...
        'LineWidth', 2.5, ...
        'Label', 'Actual', ...
        'LabelVerticalAlignment', 'bottom');

    xline(medianHneg50(validInd), 'k--', ...
        'LineWidth', 1.5, ...
        'Label', 'Null median');

    xlim([-1 1]);
    xlabel('Lag imbalance');
    ylabel('Proportion of permutations');

    title(sprintf( ...
        'H-50 ms\np = %.4f | H0/H-50 = %.3f', ...
        pHneg50(validInd), ...
        ratioH0Hneg50(validInd)));

    grid on;
    box off;

    %% ---------- overall title ----------

    if isfield(R, 'nActualSigIntPyr')
        targetPairs = R.nActualSigIntPyr;
    else
        targetPairs = NaN;
    end

    sgtitle(sprintf( ...
        ['%s | actual imbalance = %.4f | ' ...
         'real significant int-pyr pairs = %g'], ...
        animalID, ...
        observed, ...
        targetPairs), ...
        'Interpreter', 'none', ...
        'FontWeight', 'bold');

    %% ---------- save ----------

    if saveFigures

        figureBaseName = fullfile( ...
            outputFolder, ...
            sprintf('%s_lagImbalanceDistributions', animalID));

        exportgraphics(fig, ...
            [figureBaseName '.png'], ...
            'Resolution', 300);

        savefig(fig, [figureBaseName '.fig']);
    end

    if closeIndividualFigures
        close(fig);
    end

end

%% ========================================================================
%% combined distribution summary
%% ========================================================================

summaryFig = figure( ...
    'Name', 'Lag imbalance summary across animals', ...
    'Color', 'w', ...
    'Position', [100 100 1300 650]);

hold on;

xPositions = 1:nValidSessions;
horizontalOffset = 0.18;

% H0
errorbar( ...
    xPositions - horizontalOffset, ...
    medianH0, ...
    medianH0 - ciH0(:,1), ...
    ciH0(:,2) - medianH0, ...
    'o', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 7, ...
    'DisplayName', 'H0');

% H+50
errorbar( ...
    xPositions, ...
    medianH50, ...
    medianH50 - ciH50(:,1), ...
    ciH50(:,2) - medianH50, ...
    'o', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 7, ...
    'DisplayName', 'H+50 ms');

% H-50
errorbar( ...
    xPositions + horizontalOffset, ...
    medianHneg50, ...
    medianHneg50 - ciHneg50(:,1), ...
    ciHneg50(:,2) - medianHneg50, ...
    'o', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 7, ...
    'DisplayName', 'H-50 ms');

% Actual observed values
plot( ...
    xPositions, ...
    actualImbalance, ...
    'rd', ...
    'MarkerSize', 9, ...
    'LineWidth', 2, ...
    'DisplayName', 'Actual');

yline(0, 'k:');

xlim([0.5 nValidSessions + 0.5]);
ylim([-1 1]);

xticks(xPositions);
xticklabels(animalIDs);

xlabel('Animal');
ylabel('Lag imbalance');

title({ ...
    'Observed and permutation lag-imbalance values', ...
    'Points = null medians; error bars = 95% null intervals'});

legend('Location', 'best');
grid on;
box off;

if saveFigures

    exportgraphics( ...
        summaryFig, ...
        fullfile(outputFolder, ...
        'allAnimals_lagImbalanceSummary.png'), ...
        'Resolution', 300);

    savefig( ...
        summaryFig, ...
        fullfile(outputFolder, ...
        'allAnimals_lagImbalanceSummary.fig'));
end

%% ========================================================================
%% evidence-ratio summary
%% ========================================================================

evidenceFig = figure( ...
    'Name', 'Evidence ratios across animals', ...
    'Color', 'w', ...
    'Position', [100 100 1100 600]);

hold on;

plot( ...
    xPositions, ...
    ratioH0H50, ...
    'o-', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 8, ...
    'DisplayName', 'H0 / H+50');

plot( ...
    xPositions, ...
    ratioH0Hneg50, ...
    's-', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 8, ...
    'DisplayName', 'H0 / H-50');

yline(1, 'k--', ...
    'LineWidth', 1.5, ...
    'Label', 'Equal evidence');

xlim([0.5 nValidSessions + 0.5]);

xticks(xPositions);
xticklabels(animalIDs);

xlabel('Animal');
ylabel('Evidence ratio');

title('Bayes-style evidence ratios');
legend('Location', 'best');
grid on;
box off;

if saveFigures

    exportgraphics( ...
        evidenceFig, ...
        fullfile(outputFolder, ...
        'allAnimals_evidenceRatios.png'), ...
        'Resolution', 300);

    savefig( ...
        evidenceFig, ...
        fullfile(outputFolder, ...
        'allAnimals_evidenceRatios.fig'));
end

%% ========================================================================
%% print numerical summary table
%% ========================================================================

summaryTable = table( ...
    animalIDs, ...
    actualImbalance, ...
    medianH0, ...
    medianH50, ...
    medianHneg50, ...
    pH0, ...
    pH50, ...
    pHneg50, ...
    ratioH0H50, ...
    ratioH0Hneg50, ...
    selectedCountMedianH0, ...
    selectedCountMedianH50, ...
    selectedCountMedianHneg50, ...
    exactMatchRateH0, ...
    exactMatchRateH50, ...
    exactMatchRateHneg50, ...
    'VariableNames', { ...
        'Animal', ...
        'ActualImbalance', ...
        'MedianH0', ...
        'MedianH50', ...
        'MedianHneg50', ...
        'pH0', ...
        'pH50', ...
        'pHneg50', ...
        'EvidenceH0OverH50', ...
        'EvidenceH0OverHneg50', ...
        'MedianSelectedPairsH0', ...
        'MedianSelectedPairsH50', ...
        'MedianSelectedPairsHneg50', ...
        'ExactMatchRateH0', ...
        'ExactMatchRateH50', ...
        'ExactMatchRateHneg50'});

fprintf('\n========================================\n');
fprintf('LAG-IMBALANCE PLOTTING SUMMARY\n');
fprintf('========================================\n');

disp(summaryTable);

summaryCsvFile = fullfile( ...
    outputFolder, ...
    'lagImbalancePlottingSummary.csv');

writetable(summaryTable, summaryCsvFile);

fprintf('\nFigures and summary table saved to:\n%s\n', outputFolder);
fprintf('========================================\n');

end
