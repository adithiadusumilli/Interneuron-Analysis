function paperPairwiseBayesFactors(noChunkBayesFile, chunkBayesFile)
% Plot saved pairwise lag-imbalance Bayes/evidence ratios (no-chunk and chunked).
% D043 is excluded. Animals w no positive finite Bayes factors are omitted separately from each figure. Both use a common logarithmic y-axis.

% RUN:
%   paperPairwiseBayesFactors

arguments
    noChunkBayesFile (1,1) string = ...
        "C:\Users\mirilab\Documents\GlobusTransfer\pairwiseNoChunkLagImbalanceBayes50ms_StoreyCorrThresh_DAVID_UPDATE.mat"
    chunkBayesFile (1,1) string = ...
        "C:\Users\mirilab\Documents\GlobusTransfer\pairwiseChunkedLagImbalanceBayes50ms_StoreyCorrThresh_DAVID_UPDATE.mat"
end

saveDir = "X:\David\AnalysesData\InterneuronAnalyses\AA Paper Plots\pairwise bayes factors";
if ~isfolder(saveDir)
    mkdir(saveDir);
end

noChunkResults = loadResultsStruct(noChunkBayesFile,"no-chunk");
chunkResults = loadResultsStruct(chunkBayesFile,"chunked");

preferredAnimalOrder = ["D026";"D020";"D024";"D050";"D054"];
[noAnimals,noH50,noHneg50] = extractEvidenceRatios(noChunkResults);
[chAnimals,chH50,chHneg50] = extractEvidenceRatios(chunkResults);

% Exclude D043 from both datasets, including any unexpected ordering.
keepNo = ~strcmpi(noAnimals,"D043");
keepCh = ~strcmpi(chAnimals,"D043");
noAnimals = noAnimals(keepNo);
noH50 = noH50(keepNo);
noHneg50 = noHneg50(keepNo);
chAnimals = chAnimals(keepCh);
chH50 = chH50(keepCh);
chHneg50 = chHneg50(keepCh);

allAnimals = unique([noAnimals(:);chAnimals(:)],'stable');
animalOrder = strings(0,1);
for i = 1:numel(preferredAnimalOrder)
    if any(strcmpi(allAnimals,preferredAnimalOrder(i)))
        animalOrder(end+1,1) = preferredAnimalOrder(i); %#ok<AGROW>
    end
end
% Keep unexpected animal IDs rather than silently dropping them.
remaining = allAnimals(~ismember(upper(allAnimals),upper(animalOrder)));
animalOrder = [animalOrder;remaining];
if isempty(animalOrder)
    error('No animals remain after excluding D043.');
end

noH50Ordered = alignValues(noAnimals,noH50,animalOrder);
noHneg50Ordered = alignValues(noAnimals,noHneg50,animalOrder);
chH50Ordered = alignValues(chAnimals,chH50,animalOrder);
chHneg50Ordered = alignValues(chAnimals,chHneg50,animalOrder);

% A logarithmic axis requires strictly positive values. Preserve raw
% values in the printed tables; only nonpositive entries are hidden.
allRatios = [noH50Ordered(:);noHneg50Ordered(:);chH50Ordered(:);chHneg50Ordered(:)];
bad = isfinite(allRatios) & allRatios <= 0;
if any(bad)
    warning('PairwiseBayes:NonpositiveValues', ...
        '%d nonpositive evidence ratios cannot be shown on a log axis; they remain in printed tables.',sum(bad));
end
positiveValues = [allRatios(isfinite(allRatios) & allRatios > 0);1];
if numel(positiveValues) == 1
    error('No finite positive evidence ratios were found in either analysis.');
end
% Multiplicative padding is appropriate for log scales.
logMin = log10(min(positiveValues));
logMax = log10(max(positiveValues));
span = logMax-logMin;
pad = max(0.12,0.08*span);
commonYLim = [1e-1, 10^(logMax+pad)];

makeEvidenceRatioPlot(animalOrder,noH50Ordered,noHneg50Ordered, ...
    commonYLim,'No-chunk pairwise lag-imbalance Bayes factors',saveDir, ...
    'PairwiseBayesFactors_NoChunk');
makeEvidenceRatioPlot(animalOrder,chH50Ordered,chHneg50Ordered, ...
    commonYLim,'Chunked pairwise lag-imbalance Bayes factors',saveDir, ...
    'PairwiseBayesFactors_Chunked');

fprintf('\n============================================================\n');
fprintf('NO-CHUNK PAIRWISE BAYES FACTORS (D043 EXCLUDED)\n');
fprintf('============================================================\n');
disp(table(animalOrder,noH50Ordered,noHneg50Ordered, ...
    'VariableNames',{'Animal','Evidence_H0_over_H50','Evidence_H0_over_Hneg50'}));
fprintf('\n============================================================\n');
fprintf('CHUNKED PAIRWISE BAYES FACTORS (D043 EXCLUDED)\n');
fprintf('============================================================\n');
disp(table(animalOrder,chH50Ordered,chHneg50Ordered, ...
    'VariableNames',{'Animal','Evidence_H0_over_H50','Evidence_H0_over_Hneg50'}));
fprintf('Vector PDFs and editable FIG files saved to:\n%s\n',saveDir);
end

function results = loadResultsStruct(fileName,analysisLabel)
if ~isfile(fileName)
    error('%s Bayes file not found:\n%s',analysisLabel,fileName);
end
S = load(fileName,'results');
if ~isfield(S,'results')
    error('%s file does not contain variable "results":\n%s',analysisLabel,fileName);
end
results = S.results;
if ~isfield(results,'sessions') || isempty(results.sessions)
    error('%s results.sessions is missing or empty.',analysisLabel);
end
end

function [animals,h50,hneg50] = extractEvidenceRatios(results)
sessions = results.sessions;
nSess = numel(sessions);
animals = strings(nSess,1);
h50 = nan(nSess,1);
hneg50 = nan(nSess,1);
for i = 1:nSess
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
    if isfield(R,'evidenceRatio_H0_over_H50') && ~isempty(R.evidenceRatio_H0_over_H50)
        h50(i) = double(R.evidenceRatio_H0_over_H50);
    end
    if isfield(R,'evidenceRatio_H0_over_Hneg50') && ~isempty(R.evidenceRatio_H0_over_Hneg50)
        hneg50(i) = double(R.evidenceRatio_H0_over_Hneg50);
    end
end
end

function orderedVals = alignValues(sourceAnimals,sourceVals,animalOrder)
orderedVals = nan(numel(animalOrder),1);
for a = 1:numel(animalOrder)
    idx = find(strcmpi(sourceAnimals,animalOrder(a)),1,'first');
    if ~isempty(idx)
        orderedVals(a) = sourceVals(idx);
    end
end
end

function makeEvidenceRatioPlot(animalOrder,h50,hneg50,commonYLim,plotTitle,saveDir,fileStem)
% Each figure has its own x-axis animals. Keep an animal if at least one
% of its two Bayes factors is finite and strictly positive (log-plottable).
keep = (isfinite(h50) & h50 > 0) | (isfinite(hneg50) & hneg50 > 0);
animalOrder = animalOrder(keep);
h50 = h50(keep);
hneg50 = hneg50(keep);
nAnimals = numel(animalOrder);
if nAnimals == 0
    warning('PairwiseBayes:NoPlottableAnimals', ...
        'No positive finite Bayes factors for %s; figure not created.',plotTitle);
    return;
end
fig = figure('Name',plotTitle,'Color','w','Position',[100 100 720 430]);
ax = axes('Parent',fig);
hold(ax,'on');
x = (1:nAnimals)';
xOffset = 0.10;

% Keep the original within-animal connectors, only for log-plottable pairs.
for a = 1:nAnimals
    if isfinite(h50(a)) && isfinite(hneg50(a)) && h50(a)>0 && hneg50(a)>0
        line(ax,[x(a)-xOffset,x(a)+xOffset],[h50(a),hneg50(a)], ...
            'Color',[0.78 0.78 0.78],'LineWidth',1,'HandleVisibility','off');
    end
end

plotH50 = h50;
plotHneg50 = hneg50;
plotH50(~isfinite(plotH50) | plotH50<=0) = NaN;
plotHneg50(~isfinite(plotHneg50) | plotHneg50<=0) = NaN;

% Both comparisons use circles; color distinguishes the shifted model.
hPlus = scatter(ax,x-xOffset,plotH50,55,[0 0.2 0.8],'filled','o', ...
    'DisplayName','H0 / H+50');
hMinus = scatter(ax,x+xOffset,plotHneg50,55,[0.95 0.45 0.35],'filled','o', ...
    'DisplayName','H0 / H-50');
hEqual = yline(ax,1,'k--','LineWidth',1.1,'DisplayName','BF = 1 (equal evidence)');

set(ax,'YScale','log','FontSize',10,'LineWidth',0.5, ...
    'TickDir','out','TickLength',[0.015 0.015], ...
    'Box','off','XColor','k','YColor','k','Color','w');
xlim(ax,[0.5,nAnimals+0.5]);
ylim(ax,commonYLim);
xticks(ax,1:nAnimals);
xticklabels(ax,cellstr(animalOrder));
xlabel(ax,'Animal');
ylabel(ax,'Bayes factor / evidence ratio (log scale)');
title(ax,plotTitle,'FontWeight','normal','FontSize',10);
legend(ax,[hPlus,hMinus,hEqual], ...
    {'H0 / H+50','H0 / H-50','BF = 1 (equal evidence)'}, ...
    'Location','best','Box','off','FontSize',10);
grid(ax,'off');
set(fig,'Color','w');

exportgraphics(fig,fullfile(saveDir,[fileStem '.pdf']), ...
    'ContentType','vector','BackgroundColor','white');
savefig(fig,fullfile(saveDir,[fileStem '.fig']));
end
