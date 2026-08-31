function plotSavedPopAvgBayesFactors(noChunkBayesFile, chunkBayesFile)
% Plot saved overall population-average Bayes/evidence-ratio outputs.
%
% Makes exactly TWO figures:
%   1) NON-CHUNKED population-average Bayes factors
%   2) CHUNKED population-average Bayes factors
%
% In each figure:
%   x-axis = animal
%   y-axis = Bayes factor / evidence ratio
%
% Each animal has TWO points:
%   circle   = H0 / H+50
%   triangle = H0 / H-50
%
% A horizontal dashed line at Bayes factor = 1 marks equal evidence:
%   > 1 favors H0
%   < 1 favors the corresponding shifted model
%
% The same y-axis limits are used for both figures so the chunked and
% non-chunked results can be compared visually.
%
% Defaults correspond to the saved outputs from:
%   computePopAvgPeakLagBayesPosthoc50ms
%   computeChunkedPopAvgPeakLagBayesPosthoc50ms
%
% RUN:
%   plotSavedPopAvgBayesFactors
%
% OR:
%   plotSavedPopAvgBayesFactors(noChunkBayesFile, chunkBayesFile)

arguments
    noChunkBayesFile (1,1) string = ...
        "X:\David\AnalysesData\InterneuronAnalyses\Lab Meeting Pres\6 animals run cross correlation, no chunking, pop-wise\runCrossCorrelation_savedOutputs_all6Animals_PeakLagBayesPosthoc50ms.mat"

    chunkBayesFile (1,1) string = ...
        "C:\Users\mirilab\Documents\GlobusTransfer\concatCrossCorr_trialavg_chunked_popavg_ALL_SESSIONS_PeakLagBayesPosthoc50ms.mat"
end

%% ---------------- load summary tables ----------------

Tno = loadBayesSummary(noChunkBayesFile, "non-chunked");
Tch = loadBayesSummary(chunkBayesFile, "chunked");

%% ---------------- canonical animal order ----------------

preferredOrder = [ ...
    "D026", ...
    "D020", ...
    "D024", ...
    "D043", ...
    "D050", ...
    "D054"];

allAnimals = unique([string(Tno.Animal); string(Tch.Animal)],'stable');

animalOrder = strings(0,1);

for i = 1:numel(preferredOrder)
    if any(allAnimals == preferredOrder(i))
        animalOrder(end+1,1) = preferredOrder(i); %#ok<AGROW>
    end
end

% Do not silently discard any unexpected animal IDs.
remaining = allAnimals(~ismember(allAnimals,animalOrder));
animalOrder = [animalOrder; remaining];

%% ---------------- arrange values in common order ----------------

[noH50, noHneg50] = getOrderedEvidenceRatios(Tno,animalOrder);
[chH50, chHneg50] = getOrderedEvidenceRatios(Tch,animalOrder);

%% ---------------- common y-axis limits ----------------

allVals = [ ...
    noH50(:); ...
    noHneg50(:); ...
    chH50(:); ...
    chHneg50(:); ...
    1];

allVals = allVals(isfinite(allVals));

if isempty(allVals)
    error('No finite Bayes factors were found in either file.');
end

yMin = min(allVals);
yMax = max(allVals);

if yMin == yMax
    yPad = max(0.1,0.1*abs(yMin));
else
    yPad = 0.08*(yMax-yMin);
end

commonYLim = [max(0,yMin-yPad), yMax+yPad];

% Ensure BF=1 is comfortably visible.
commonYLim(1) = min(commonYLim(1),0.9);
commonYLim(2) = max(commonYLim(2),1.1);

%% ---------------- figure 1: non-chunk ----------------

makeBayesDotPlot( ...
    animalOrder, ...
    noH50, ...
    noHneg50, ...
    commonYLim, ...
    'Non-chunked population-average Bayes factors');

%% ---------------- figure 2: chunk ----------------

makeBayesDotPlot( ...
    animalOrder, ...
    chH50, ...
    chHneg50, ...
    commonYLim, ...
    'Chunked trial-averaged population Bayes factors');

%% ---------------- print values ----------------

fprintf('\n============================================================\n');
fprintf('NON-CHUNKED POPULATION-AVERAGE BAYES FACTORS\n');
fprintf('============================================================\n');

disp(table( ...
    animalOrder, ...
    noH50, ...
    noHneg50, ...
    'VariableNames',{ ...
        'Animal', ...
        'Evidence_H0_over_H50', ...
        'Evidence_H0_over_Hneg50'}));

fprintf('\n============================================================\n');
fprintf('CHUNKED POPULATION-AVERAGE BAYES FACTORS\n');
fprintf('============================================================\n');

disp(table( ...
    animalOrder, ...
    chH50, ...
    chHneg50, ...
    'VariableNames',{ ...
        'Animal', ...
        'Evidence_H0_over_H50', ...
        'Evidence_H0_over_Hneg50'}));

end

%% ========================================================================
%% helper: load Bayes summary table
%% ========================================================================

function T = loadBayesSummary(fileName,analysisLabel)

if ~isfile(fileName)
    error('%s Bayes output file not found:\n%s',analysisLabel,fileName);
end

S = load(fileName,'summaryTable','posthocResults');

if isfield(S,'summaryTable') && ~isempty(S.summaryTable)
    T = S.summaryTable;

elseif isfield(S,'posthocResults') && ...
        isfield(S.posthocResults,'summaryTable') && ...
        ~isempty(S.posthocResults.summaryTable)

    T = S.posthocResults.summaryTable;

else
    error('%s Bayes file does not contain a summaryTable.',analysisLabel);
end

requiredVars = { ...
    'Animal', ...
    'Evidence_H0_over_H50', ...
    'Evidence_H0_over_Hneg50'};

for i = 1:numel(requiredVars)
    if ~ismember(requiredVars{i},T.Properties.VariableNames)
        error('%s summaryTable is missing "%s".', ...
            analysisLabel,requiredVars{i});
    end
end

T.Animal = string(T.Animal);

end

%% ========================================================================
%% helper: align table values to requested animal order
%% ========================================================================

function [h50,hneg50] = getOrderedEvidenceRatios(T,animalOrder)

nAnimals = numel(animalOrder);

h50 = nan(nAnimals,1);
hneg50 = nan(nAnimals,1);

for a = 1:nAnimals

    idx = find(T.Animal == animalOrder(a),1,'first');

    if isempty(idx)
        continue;
    end

    h50(a) = double(T.Evidence_H0_over_H50(idx));
    hneg50(a) = double(T.Evidence_H0_over_Hneg50(idx));
end

end

%% ========================================================================
%% helper: one dot plot
%% ========================================================================

function makeBayesDotPlot( ...
    animalOrder, ...
    evidenceH50, ...
    evidenceHneg50, ...
    commonYLim, ...
    plotTitle)

nAnimals = numel(animalOrder);

figure( ...
    'Name',plotTitle, ...
    'Color','w');

hold on;

x = (1:nAnimals)';
offset = 0.10;

% Light within-animal connector between the two model comparisons.
for a = 1:nAnimals

    if isfinite(evidenceH50(a)) && isfinite(evidenceHneg50(a))

        line( ...
            [x(a)-offset, x(a)+offset], ...
            [evidenceH50(a), evidenceHneg50(a)], ...
            'Color',[0.78 0.78 0.78], ...
            'LineWidth',1);
    end
end

hPlus = scatter( ...
    x-offset, ...
    evidenceH50, ...
    85, ...
    'filled', ...
    'Marker','o');

hMinus = scatter( ...
    x+offset, ...
    evidenceHneg50, ...
    85, ...
    'filled', ...
    'Marker','^');

yline(1,'k--','LineWidth',1.5);

xlim([0.5 nAnimals+0.5]);
ylim(commonYLim);

xticks(1:nAnimals);
xticklabels(cellstr(animalOrder));

xlabel('Animal');
ylabel('Bayes factor / evidence ratio');

title(plotTitle);

legend( ...
    [hPlus,hMinus], ...
    {'H0 / H+50','H0 / H-50'}, ...
    'Location','best', ...
    'Box','off');

grid on;
box off;

end
