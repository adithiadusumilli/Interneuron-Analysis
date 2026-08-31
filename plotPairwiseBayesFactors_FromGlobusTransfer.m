function plotPairwiseBayesFactors_FromGlobusTransfer(noChunkBayesFile, chunkBayesFile)
% Plot saved pairwise lag-imbalance Bayes/evidence-ratio outputs OUTSIDE QUEST.
%
% This function expects the two QUEST-generated .mat files to have been
% copied into the Windows GlobusTransfer folder on C:.
%
% Makes exactly TWO figures:
%   1) NO-CHUNK pairwise lag-imbalance Bayes factors
%   2) CHUNKED pairwise lag-imbalance Bayes factors
%
% In each figure:
%   x-axis = animal
%   y-axis = Bayes factor / evidence ratio
%
% Each animal has two slightly offset dots:
%   circle   = H0 / H+50
%   triangle = H0 / H-50
%
% A horizontal line at BF = 1 marks equal evidence:
%   BF > 1 favors H0
%   BF < 1 favors the corresponding shifted model
%
% DEFAULT INPUT FILES:
%
%   C:\Users\mirilab\Documents\GlobusTransfer\
%       pairwiseNoChunkLagImbalanceBayes50ms_StoreyCorrThresh_DAVID_UPDATE.mat
%
%   C:\Users\mirilab\Documents\GlobusTransfer\
%       pairwiseChunkedLagImbalanceBayes50ms_StoreyCorrThresh_DAVID_UPDATE.mat
%
% RUN:
%   plotPairwiseBayesFactors_FromGlobusTransfer
%
% Or:
%   plotPairwiseBayesFactors_FromGlobusTransfer(noChunkBayesFile,chunkBayesFile)

arguments
    noChunkBayesFile (1,1) string = ...
        "C:\Users\mirilab\Documents\GlobusTransfer\pairwiseNoChunkLagImbalanceBayes50ms_StoreyCorrThresh_DAVID_UPDATE.mat"

    chunkBayesFile (1,1) string = ...
        "C:\Users\mirilab\Documents\GlobusTransfer\pairwiseChunkedLagImbalanceBayes50ms_StoreyCorrThresh_DAVID_UPDATE.mat"
end

%% ========================================================================
%% LOAD SAVED QUEST OUTPUTS
%% ========================================================================

noChunkResults = loadResultsStruct(noChunkBayesFile,"no-chunk");
chunkResults   = loadResultsStruct(chunkBayesFile,"chunked");

%% ========================================================================
%% CANONICAL ANIMAL ORDER
%% ========================================================================

preferredAnimalOrder = [ ...
    "D026", ...
    "D020", ...
    "D024", ...
    "D043", ...
    "D050", ...
    "D054"];

[noAnimals,noH50,noHneg50] = extractEvidenceRatios(noChunkResults);
[chAnimals,chH50,chHneg50] = extractEvidenceRatios(chunkResults);

allAnimals = unique([noAnimals(:);chAnimals(:)],'stable');

animalOrder = strings(0,1);

for i = 1:numel(preferredAnimalOrder)
    if any(allAnimals == preferredAnimalOrder(i))
        animalOrder(end+1,1) = preferredAnimalOrder(i); %#ok<AGROW>
    end
end

% Preserve unexpected IDs rather than dropping them.
remaining = allAnimals(~ismember(allAnimals,animalOrder));
animalOrder = [animalOrder;remaining];

%% ========================================================================
%% ALIGN BOTH ANALYSES TO SAME ANIMAL ORDER
%% ========================================================================

noH50Ordered    = alignValues(noAnimals,noH50,animalOrder);
noHneg50Ordered = alignValues(noAnimals,noHneg50,animalOrder);

chH50Ordered    = alignValues(chAnimals,chH50,animalOrder);
chHneg50Ordered = alignValues(chAnimals,chHneg50,animalOrder);

%% ========================================================================
%% USE COMMON Y LIMITS FOR DIRECT VISUAL COMPARISON
%% ========================================================================

allVals = [ ...
    noH50Ordered(:); ...
    noHneg50Ordered(:); ...
    chH50Ordered(:); ...
    chHneg50Ordered(:); ...
    1];

allVals = allVals(isfinite(allVals));

if isempty(allVals)
    error('No finite Bayes/evidence-ratio values were found.');
end

yMin = min(allVals);
yMax = max(allVals);

if yMin == yMax
    yPad = max(0.1,0.1*abs(yMin));
else
    yPad = 0.08*(yMax-yMin);
end

commonYLim = [max(0,yMin-yPad), yMax+yPad];

% Guarantee BF=1 is visible.
commonYLim(1) = min(commonYLim(1),0.9);
commonYLim(2) = max(commonYLim(2),1.1);

%% ========================================================================
%% FIGURE 1: NO-CHUNK
%% ========================================================================

makeEvidenceRatioPlot( ...
    animalOrder, ...
    noH50Ordered, ...
    noHneg50Ordered, ...
    commonYLim, ...
    'No-chunk pairwise lag-imbalance Bayes factors');

%% ========================================================================
%% FIGURE 2: CHUNKED
%% ========================================================================

makeEvidenceRatioPlot( ...
    animalOrder, ...
    chH50Ordered, ...
    chHneg50Ordered, ...
    commonYLim, ...
    'Chunked pairwise lag-imbalance Bayes factors');

%% ========================================================================
%% PRINT EXACT VALUES
%% ========================================================================

fprintf('\n============================================================\n');
fprintf('NO-CHUNK PAIRWISE BAYES FACTORS\n');
fprintf('============================================================\n');

noChunkTable = table( ...
    animalOrder, ...
    noH50Ordered, ...
    noHneg50Ordered, ...
    'VariableNames',{ ...
        'Animal', ...
        'Evidence_H0_over_H50', ...
        'Evidence_H0_over_Hneg50'});

disp(noChunkTable);

fprintf('\n============================================================\n');
fprintf('CHUNKED PAIRWISE BAYES FACTORS\n');
fprintf('============================================================\n');

chunkTable = table( ...
    animalOrder, ...
    chH50Ordered, ...
    chHneg50Ordered, ...
    'VariableNames',{ ...
        'Animal', ...
        'Evidence_H0_over_H50', ...
        'Evidence_H0_over_Hneg50'});

disp(chunkTable);

end

%% ========================================================================
%% HELPER: LOAD RESULTS STRUCT
%% ========================================================================

function results = loadResultsStruct(fileName,analysisLabel)

if ~isfile(fileName)
    error('%s Bayes file not found:\n%s',analysisLabel,fileName);
end

S = load(fileName,'results');

if ~isfield(S,'results')
    error('%s file does not contain variable "results":\n%s', ...
        analysisLabel,fileName);
end

results = S.results;

if ~isfield(results,'sessions') || isempty(results.sessions)
    error('%s results.sessions is missing or empty.',analysisLabel);
end

end

%% ========================================================================
%% HELPER: EXTRACT SAVED ANIMAL IDS + BAYES FACTORS
%% ========================================================================

function [animals,h50,hneg50] = extractEvidenceRatios(results)

sessions = results.sessions;
nSess = numel(sessions);

animals = strings(nSess,1);
h50 = nan(nSess,1);
hneg50 = nan(nSess,1);

for i = 1:nSess

    % Both pairwise files save results.sessions as cells.
    if iscell(sessions)
        R = sessions{i};
    else
        R = sessions(i);
    end

    if isempty(R)
        animals(i) = "Animal_" + i;
        continue;
    end

    if isfield(R,'animalID') && ~isempty(R.animalID)
        animals(i) = string(R.animalID);
    elseif isfield(R,'mouseID') && ~isempty(R.mouseID)
        animals(i) = string(R.mouseID);
    elseif isfield(R,'baseDir') && ~isempty(R.baseDir)
        hit = regexp(char(R.baseDir),'D\d+','match','once');
        if ~isempty(hit)
            animals(i) = string(hit);
        else
            animals(i) = "Animal_" + i;
        end
    else
        animals(i) = "Animal_" + i;
    end

    if isfield(R,'evidenceRatio_H0_over_H50') && ...
            ~isempty(R.evidenceRatio_H0_over_H50)

        h50(i) = double(R.evidenceRatio_H0_over_H50);
    end

    if isfield(R,'evidenceRatio_H0_over_Hneg50') && ...
            ~isempty(R.evidenceRatio_H0_over_Hneg50)

        hneg50(i) = double(R.evidenceRatio_H0_over_Hneg50);
    end
end

end

%% ========================================================================
%% HELPER: ALIGN VALUES TO REQUESTED ANIMAL ORDER
%% ========================================================================

function orderedVals = alignValues(sourceAnimals,sourceVals,animalOrder)

orderedVals = nan(numel(animalOrder),1);

for a = 1:numel(animalOrder)

    idx = find(sourceAnimals == animalOrder(a),1,'first');

    if ~isempty(idx)
        orderedVals(a) = sourceVals(idx);
    end
end

end

%% ========================================================================
%% HELPER: MAKE ONE BAYES-FACTOR DOT PLOT
%% ========================================================================

function makeEvidenceRatioPlot( ...
    animalOrder,h50,hneg50,commonYLim,plotTitle)

nAnimals = numel(animalOrder);

figure( ...
    'Name',plotTitle, ...
    'Color','w');

hold on;

x = (1:nAnimals)';
xOffset = 0.10;

% Connect the two model comparisons for the same animal.
for a = 1:nAnimals

    if isfinite(h50(a)) && isfinite(hneg50(a))

        line( ...
            [x(a)-xOffset,x(a)+xOffset], ...
            [h50(a),hneg50(a)], ...
            'Color',[0.78 0.78 0.78], ...
            'LineWidth',1);
    end
end

hPlus = scatter( ...
    x-xOffset, ...
    h50, ...
    85, ...
    'filled', ...
    'Marker','o');

hMinus = scatter( ...
    x+xOffset, ...
    hneg50, ...
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
