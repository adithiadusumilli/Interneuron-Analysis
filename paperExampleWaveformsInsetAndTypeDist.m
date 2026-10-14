function paperExampleWaveformsInsetAndTypeDist(baseDir)
% paper waveform insets and M1 neuron-type percentages
% includes D026, D020, D024, D050, D054 (excludes D043)
% uses classifications and file paths saved by classifyCortexAndStriatum_getMouseDataNames

% Example: plotExampleWaveformsForInset("Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData")

% The example animal controls the two waveform insets. The percentage
% plot always includes all five specified animals.

arguments
    baseDir (1,1) string
end

%% Settings
classificationFile = "X:\David\AnalysesData\AA_classifications.mat";
saveDir = "X:\David\AnalysesData\InterneuronAnalyses\AA Paper Plots\classifications";
if ~isfolder(saveDir), mkdir(saveDir); end
fontSize = 10;
axesLineWidth = 0.5;
tickLength = [0.025 0.025];
fs = 30000;
wanted = ["D026", "D020", "D024", "D050", "D054"];
markers = {'o','s','d','^','v'};

%% Load classifications and verify animal order
if ~isfile(classificationFile)
    error('Cannot find classification file: %s', classificationFile);
end
C = load(classificationFile);
if ~isfield(C,'classifications') || ~isfield(C,'mouseIDs')
    error('AA_classifications.mat must contain classifications and mouseIDs.');
end
savedIDs = string(C.mouseIDs(:));
rows = zeros(numel(wanted),1);
for a = 1:numel(wanted)
    hits = find(strcmpi(savedIDs,wanted(a)));
    if numel(hits) ~= 1
        error('Expected one classification entry for %s; found %d.',wanted(a),numel(hits));
    end
    rows(a) = hits;
end

match = regexp(char(baseDir),'D\d{3}','match','once');
if isempty(match)
    error('Could not identify mouse ID in baseDir: %s',baseDir);
end
exampleID = string(match);
if exampleID == "D043"
    error('D043 is excluded from this paper figure.');
end
exampleAnimal = find(wanted == exampleID,1);
if isempty(exampleAnimal)
    error('Unsupported example animal %s.',exampleID);
end

%% Percentages for five animals from saved cortex classifications
nAnimals = numel(wanted);
pyrPctAll = nan(nAnimals,1);
intPctAll = nan(nAnimals,1);
nPyrAll = zeros(nAnimals,1);
nIntAll = zeros(nAnimals,1);
nTotAll = zeros(nAnimals,1);

for a = 1:nAnimals
    labels = double(C.classifications{rows(a),1});
    labels = labels(:);
    % The saved labels are full-neuron vectors: NaN marks other regions
    % and neurons in the classification gap.
    nPyrAll(a) = sum(labels == 0);
    nIntAll(a) = sum(labels == 1);
    nTotAll(a) = nPyrAll(a) + nIntAll(a);
    if nTotAll(a) == 0
        error('No classified cortex units for %s.',wanted(a));
    end
    pyrPctAll(a) = 100*nPyrAll(a)/nTotAll(a);
    intPctAll(a) = 100*nIntAll(a)/nTotAll(a);
end

fprintf('\nM1 neuron-type percentages (classified cortex units only):\n');
for a = 1:nAnimals
    fprintf('%s: pyramidal %.2f%% (%d/%d), interneuron %.2f%% (%d/%d)\n', ...
        wanted(a),pyrPctAll(a),nPyrAll(a),nTotAll(a), ...
        intPctAll(a),nIntAll(a),nTotAll(a));
end
fprintf('Mean across animals: pyramidal %.2f%%, interneuron %.2f%%\n', ...
    mean(pyrPctAll),mean(intPctAll));

%% Load example animal waveform data, using saved David-helper path
row = rows(exampleAnimal);
exampleNeuronFile = fullfile(baseDir,'neuronDataStruct.mat');
if isfield(C,'neuronDataStructFiles') && numel(C.neuronDataStructFiles) >= row
    candidate = string(C.neuronDataStructFiles{row});
    if isfile(candidate)
        exampleNeuronFile = candidate;
    end
end
if ~isfile(exampleNeuronFile)
    error('Cannot find example neuronDataStruct: %s',exampleNeuronFile);
end
S = load(exampleNeuronFile,'neuronDataStruct');
if ~isfield(S,'neuronDataStruct')
    error('Missing neuronDataStruct in %s.',exampleNeuronFile);
end
neurons = S.neuronDataStruct;
labels = double(C.classifications{row,1});
labels = labels(:);
if numel(labels) ~= numel(neurons)
    error('Saved cortex labels (%d) do not match neuronDataStruct (%d) for %s.', ...
        numel(labels),numel(neurons),exampleID);
end
pyrInds = find(labels == 0);
intInds = find(labels == 1);
if isempty(pyrInds) || isempty(intInds)
    error('%s does not have both classified cortex neuron types.',exampleID);
end

% Use trough-to-SUBSEQUENT-peak widths, matching the classification code.
pyrExampleInd = chooseExample(neurons,pyrInds,fs,75);
intExampleInd = chooseExample(neurons,intInds,fs,25);

figPyr = plotSingleWaveformInset(neurons(pyrExampleInd),fs,[0 0 1], ...
    sprintf('%s Pyramidal Example',exampleID),fontSize,axesLineWidth,tickLength);
savePaper(figPyr,saveDir,sprintf('%s_pyramidal_waveform_inset',exampleID));

figInt = plotSingleWaveformInset(neurons(intExampleInd),fs,[1 0 0], ...
    sprintf('%s Interneuron Example',exampleID),fontSize,axesLineWidth,tickLength);
savePaper(figInt,saveDir,sprintf('%s_interneuron_waveform_inset',exampleID));

%% Combined five-animal percentage plot
figPct = figure('Name','M1 Neuron-Type Percentages','Color','w', ...
    'Position',[100 100 580 390]);
ax = axes(figPct); hold(ax,'on');
b1 = bar(ax,1,mean(pyrPctAll),'FaceColor',[0 0 1], ...
    'FaceAlpha',0.28,'EdgeColor','none','BarWidth',0.7);
b2 = bar(ax,2,mean(intPctAll),'FaceColor',[1 0 0], ...
    'FaceAlpha',0.28,'EdgeColor','none','BarWidth',0.7);

jitter = linspace(-0.07,0.07,nAnimals)';
xPyr = 1+jitter;
xInt = 2+jitter;
legendHandles = gobjects(nAnimals,1);
for a = 1:nAnimals
    plot(ax,[xPyr(a),xInt(a)],[pyrPctAll(a),intPctAll(a)],'-', ...
        'Color',[0.55 0.55 0.55],'LineWidth',0.9,'HandleVisibility','off');
    plot(ax,xPyr(a),pyrPctAll(a),markers{a},'MarkerFaceColor',[0 0 1], ...
        'MarkerEdgeColor','k','MarkerSize',5,'HandleVisibility','off');
    legendHandles(a) = plot(ax,xInt(a),intPctAll(a),markers{a}, ...
        'MarkerFaceColor',[1 0 0],'MarkerEdgeColor','k','MarkerSize',5, ...
        'LineStyle','none','DisplayName',char(wanted(a)));
end
xlim(ax,[0.4 2.6]); ylim(ax,[0 100]);
xticks(ax,[1 2]); xticklabels(ax,{'Pyramidal Neuron','Interneuron'});
ylabel(ax,'Percent of M1 Units (%)');
title(ax,'M1 Neuron-Type Percentages','FontWeight','normal');
legend(ax,legendHandles,cellstr(wanted),'Location','eastoutside', ...
    'FontSize',fontSize,'Box','off');
set(ax,'FontSize',fontSize,'LineWidth',axesLineWidth,'TickDir','out', ...
    'TickLength',tickLength,'XColor','k','YColor','k','Color','w');
grid(ax,'off'); box(ax,'off');
savePaper(figPct,saveDir,'M1_neuron_type_percentages_D026_D020_D024_D050_D054');
fprintf('All figures saved to: %s\n',saveDir);
end

function chosen = chooseExample(neurons,inds,fs,percentile)
widths = nan(numel(inds),1);
for k = 1:numel(inds)
    try
        unit = neurons(inds(k));
        ap = double(unit.waveforms(:,unit.biggestChan));
        [~,mn] = min(ap);
        [~,mxAfter] = max(ap(mn:end));
        widths(k) = (mxAfter-1)/fs;
    catch
        widths(k) = NaN;
    end
end
valid = isfinite(widths);
if ~any(valid)
    error('No valid waveform widths found for example selection.');
end
validInds = inds(valid);
validWidths = widths(valid);
target = prctile(validWidths,percentile);
[~,i] = min(abs(validWidths-target));
chosen = validInds(i);
end

function fig = plotSingleWaveformInset(unit,fs,waveColor,figTitle,fontSize,axesLineWidth,tickLength)
ap = double(unit.waveforms(:,unit.biggestChan));
[~,mn] = min(ap);
[~,mxAfter] = max(ap(mn:end));
mx = mn+mxAfter-1;
if mx <= mn
    error('Example waveform has no subsequent peak after trough.');
end
leftPad = 8;
rightPad = 12;
i0 = max(1,mn-leftPad);
i1 = min(numel(ap),mx+rightPad);
inds = i0:i1;
tMs = (inds-1)/fs*1000;
y = ap(inds);
fig = figure('Name',figTitle,'Color','w','Position',[100 100 310 240]);
ax = axes(fig); hold(ax,'on');
plot(ax,tMs,y,'Color',waveColor,'LineWidth',1.5);
plot(ax,(mn-1)/fs*1000,ap(mn),'o','MarkerSize',4, ...
    'MarkerFaceColor',waveColor,'MarkerEdgeColor',waveColor);
plot(ax,(mx-1)/fs*1000,ap(mx),'o','MarkerSize',4, ...
    'MarkerFaceColor',waveColor,'MarkerEdgeColor',waveColor);
xlim(ax,[tMs(1),tMs(end)]);
yRange = max(y)-min(y);
if yRange <= 0, yRange = max(abs(y))*0.2+eps; end
ylim(ax,[min(y)-0.25*yRange,max(y)+0.15*yRange]);
% Only draw a true 0.5 ms bar when it fits in the displayed x-range.
xRange = tMs(end)-tMs(1);
if xRange >= 0.65
    x0 = tMs(1)+0.08*xRange;
    yLim = ylim(ax);
    y0 = yLim(1)+0.18*diff(yLim);
    plot(ax,[x0,x0+0.5],[y0,y0],'k','LineWidth',1.2);
    text(ax,x0+0.25,y0-0.07*diff(yLim),'0.5 ms', ...
        'HorizontalAlignment','center','VerticalAlignment','top','FontSize',fontSize);
end
title(ax,figTitle,'FontSize',fontSize,'FontWeight','normal');
set(ax,'FontSize',fontSize,'LineWidth',axesLineWidth,'TickDir','out', ...
    'TickLength',tickLength,'Color','w','XColor','none','YColor','none');
grid(ax,'off'); box(ax,'off');
end

function savePaper(fig,saveDir,stem)
pdfPath = fullfile(saveDir,[stem '.pdf']);
figPath = fullfile(saveDir,[stem '.fig']);
exportgraphics(fig,pdfPath,'ContentType','vector','BackgroundColor','white');
savefig(fig,figPath);
fprintf('Saved: %s\n',pdfPath);
end
