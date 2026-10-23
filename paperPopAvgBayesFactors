function paperPopAvgBayesFactors(noChunkBayesFile, chunkBayesFile)
% Plot saved population-average Bayes/evidence ratios for five animals.
% Makes two figures: non-chunked and chunked. D043 is excluded.
% Both +/-50 ms comparisons use circles distinguished by color.
% A dashed BF=1 line marks equal evidence; axes use a shared log scale.
%
% RUN:
%   paperPopAvgBayesFactors
%   paperPopAvgBayesFactors(noChunkBayesFile, chunkBayesFile)

arguments
    noChunkBayesFile (1,1) string = ...
        "X:\David\AnalysesData\InterneuronAnalyses\Lab Meeting Pres\6 animals run cross correlation, no chunking, pop-wise\runCrossCorrelation_savedOutputs_all6Animals_PeakLagBayesPosthoc50ms.mat"
    chunkBayesFile (1,1) string = ...
        "C:\Users\mirilab\Documents\GlobusTransfer\concatCrossCorr_trialavg_chunked_popavg_ALL_SESSIONS_PeakLagBayesPosthoc50ms.mat"
end

saveDir = "X:\David\AnalysesData\InterneuronAnalyses\AA Paper Plots\popavg bayes factors";
if ~isfolder(saveDir)
    mkdir(saveDir);
end

%% Load saved Bayes summary tables
Tno = loadBayesSummary(noChunkBayesFile, "non-chunked");
Tch = loadBayesSummary(chunkBayesFile, "chunked");

%% Canonical five-animal order (exclude D043)
preferredOrder = ["D026"; "D020"; "D024"; "D050"; "D054"];
allAnimals = unique([string(Tno.Animal); string(Tch.Animal)], 'stable');
allAnimals = allAnimals(allAnimals ~= "D043");

animalOrder = preferredOrder(ismember(preferredOrder, allAnimals));
% Retain unexpected animals (except D043), as in the original function.
remaining = allAnimals(~ismember(allAnimals, animalOrder));
animalOrder = [animalOrder; remaining];

if isempty(animalOrder)
    error('No animals other than D043 were found in the saved Bayes summaries.');
end

%% Align ratios by animal
[noH50, noHneg50] = getOrderedEvidenceRatios(Tno, animalOrder);
[chH50, chHneg50] = getOrderedEvidenceRatios(Tch, animalOrder);

%% Shared positive log-axis limits for comparison
allVals = [noH50(:); noHneg50(:); chH50(:); chHneg50(:); 1];
invalid = isfinite(allVals) & allVals <= 0;
if any(invalid)
    warning('plotSavedPopAvgBayesFactors:NonpositiveRatios', ...
        '%d finite Bayes ratios are <= 0 and cannot appear on a log axis.', sum(invalid));
end
validVals = allVals(isfinite(allVals) & allVals > 0);
if isempty(validVals)
    error('No positive finite Bayes factors were found.');
end
% Pad in log space (multiplicative rather than additive padding).
logMin = log10(min(validVals));
logMax = log10(max(validVals));
logRange = logMax - logMin;
logPad = max(0.10, 0.08 * logRange);
commonYLim = 10.^[logMin - logPad, logMax + logPad];

%% Figure 1: non-chunked
makeBayesDotPlot(animalOrder, noH50, noHneg50, commonYLim, ...
    'Non-chunked population-average Bayes factors', saveDir, ...
    'PopAvgBayesFactors_NonChunked_D043Excluded');

%% Figure 2: chunked
makeBayesDotPlot(animalOrder, chH50, chHneg50, commonYLim, ...
    'Chunked trial-averaged population Bayes factors', saveDir, ...
    'PopAvgBayesFactors_Chunked_D043Excluded');

%% Print the same values as before, excluding D043
fprintf('\n============================================================\n');
fprintf('NON-CHUNKED POPULATION-AVERAGE BAYES FACTORS (D043 EXCLUDED)\n');
fprintf('============================================================\n');
disp(table(animalOrder, noH50, noHneg50, 'VariableNames', ...
    {'Animal','Evidence_H0_over_H50','Evidence_H0_over_Hneg50'}));

fprintf('\n============================================================\n');
fprintf('CHUNKED POPULATION-AVERAGE BAYES FACTORS (D043 EXCLUDED)\n');
fprintf('============================================================\n');
disp(table(animalOrder, chH50, chHneg50, 'VariableNames', ...
    {'Animal','Evidence_H0_over_H50','Evidence_H0_over_Hneg50'}));
fprintf('Saved vector PDFs and MATLAB FIG files to: %s\n', saveDir);
end

%% Helper: load Bayes summary table
function T = loadBayesSummary(fileName, analysisLabel)
if ~isfile(fileName)
    error('%s Bayes output file not found:\n%s', analysisLabel, fileName);
end
S = load(fileName, 'summaryTable', 'posthocResults');
if isfield(S, 'summaryTable') && ~isempty(S.summaryTable)
    T = S.summaryTable;
elseif isfield(S, 'posthocResults') && ...
        isfield(S.posthocResults, 'summaryTable') && ...
        ~isempty(S.posthocResults.summaryTable)
    T = S.posthocResults.summaryTable;
else
    error('%s Bayes file does not contain a summaryTable.', analysisLabel);
end
requiredVars = {'Animal','Evidence_H0_over_H50','Evidence_H0_over_Hneg50'};
for i = 1:numel(requiredVars)
    if ~ismember(requiredVars{i}, T.Properties.VariableNames)
        error('%s summaryTable is missing "%s".', analysisLabel, requiredVars{i});
    end
end
T.Animal = string(T.Animal);
end

%% Helper: align table values to animal order
function [h50, hneg50] = getOrderedEvidenceRatios(T, animalOrder)
nAnimals = numel(animalOrder);
h50 = nan(nAnimals,1);
hneg50 = nan(nAnimals,1);
for a = 1:nAnimals
    idx = find(T.Animal == animalOrder(a), 1, 'first');
    if isempty(idx), continue; end
    h50(a) = double(T.Evidence_H0_over_H50(idx));
    hneg50(a) = double(T.Evidence_H0_over_Hneg50(idx));
end
end

%% Helper: paper-style dot plot
function makeBayesDotPlot(animalOrder, evidenceH50, evidenceHneg50, ...
    commonYLim, plotTitle, saveDir, fileStem)

nAnimals = numel(animalOrder);
fig = figure('Name', plotTitle, 'Color', 'w', ...
    'Position', [100 100 680 410]);
ax = axes('Parent', fig);
hold(ax, 'on');

x = (1:nAnimals)';
offset = 0.10;

% Logarithmic axes only permit positive finite ratios.
plus = evidenceH50;
minus = evidenceHneg50;
plus(~isfinite(plus) | plus <= 0) = NaN;
minus(~isfinite(minus) | minus <= 0) = NaN;

for a = 1:nAnimals
    if isfinite(plus(a)) && isfinite(minus(a))
        line(ax, [x(a)-offset, x(a)+offset], [plus(a), minus(a)], ...
            'Color', [0.78 0.78 0.78], 'LineWidth', 0.8, ...
            'HandleVisibility', 'off');
    end
end

% Same circular marker for both comparisons; distinguish by color.
hPlus = scatter(ax, x-offset, plus, 42, [0 0.2 0.6], ...
    'o', 'filled', 'DisplayName', 'H0 / H+50');
hMinus = scatter(ax, x+offset, minus, 42, [0.95 0.45 0.35], ...
    'o', 'filled', 'DisplayName', 'H0 / H-50');
hEqual = yline(ax, 1, 'k--', 'LineWidth', 1, ...
    'DisplayName', 'BF = 1 (equal evidence)');

set(ax, 'YScale', 'log', 'FontSize', 10, 'LineWidth', 0.5, ...
    'TickDir', 'out', 'TickLength', [0.015 0.015], ...
    'XColor', 'k', 'YColor', 'k', 'Color', 'w', 'Box', 'off');
xlim(ax, [0.5 nAnimals+0.5]);
ylim(ax, commonYLim);
xticks(ax, 1:nAnimals);
xticklabels(ax, cellstr(animalOrder));
xlabel(ax, 'Animal');
ylabel(ax, 'Bayes factor / evidence ratio (log scale)');
title(ax, plotTitle, 'FontSize', 10, 'FontWeight', 'normal');
legend(ax, [hPlus hMinus hEqual], ...
    {'H0 / H+50','H0 / H-50','BF = 1 (equal evidence)'}, ...
    'Location', 'best', 'Box', 'off', 'FontSize', 10);
grid(ax, 'off');

exportgraphics(fig, fullfile(saveDir, fileStem + ".pdf"), ...
    'ContentType', 'vector', 'BackgroundColor', 'white');
savefig(fig, fullfile(saveDir, fileStem + ".fig"));
end
