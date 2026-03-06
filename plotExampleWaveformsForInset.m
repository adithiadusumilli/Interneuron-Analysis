function plotExampleWaveformsForInset(baseDir)
% plots one example cortex pyramidal waveform and one example cortex interneuron waveform
% in separate small figures for use as insets

% j run:
% plotExampleWaveformsForInset("Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData")

arguments
    baseDir (1,1) string
end

%% ---- settings ----
conslidatedDataFoler = 'X:\David\AnalysesData';

animalFolders = {
    'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData', ...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData', ...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData', ...
    'X:\David\ArenaRecordings\D043-020425-ArenaRecording\ProcessedData'
};

fs = 30000; % hz

%% ---- load data ----
S = load(fullfile(baseDir,'neuronDataStruct.mat'),'neuronDataStruct');
F = load(fullfile(baseDir,'NeuralFiringRates1msBins10msGauss.mat'),'cortexInds');
C = load(fullfile(conslidatedDataFoler,'AA_classifications.mat'),'classifications');

neuronDataStruct = S.neuronDataStruct;
cortexInds = double(F.cortexInds(:))';
classifications = C.classifications;

matchRow = find(contains(string(animalFolders), baseDir), 1);
if isempty(matchRow)
    error('could not match baseDir to animalFolders.');
end

cortexLabelsAll = double(classifications{matchRow,1});
cortexLabels = cortexLabelsAll(cortexInds); % 0 = pyr, 1 = int

%% ---- pick one example pyr and one example int ----
pyrInds = cortexInds(cortexLabels == 0);
intInds = cortexInds(cortexLabels == 1);

if isempty(pyrInds) || isempty(intInds)
    error('could not find both pyramidal and interneuron cortex units.');
end

% choose median-width examples rather than arbitrary first units
[pyrExampleInd, intExampleInd] = pickRepresentativeExamples(neuronDataStruct, pyrInds, intInds, fs);

%% ---- plot pyramidal waveform ----
plotSingleWaveformInset(neuronDataStruct(pyrExampleInd), fs, [0 0 1], 'Pyramidal Example');

%% ---- plot interneuron waveform ----
plotSingleWaveformInset(neuronDataStruct(intExampleInd), fs, [1 0 0], 'Interneuron Example');

end

function [pyrExampleInd, intExampleInd] = pickRepresentativeExamples(neuronDataStruct, pyrInds, intInds, fs)
% choose a representative unit near the median width for each class

pyrWidths = nan(1,numel(pyrInds));
for i = 1:numel(pyrInds)
    waveform = neuronDataStruct(pyrInds(i)).waveforms;
    biggestChan = neuronDataStruct(pyrInds(i)).biggestChan;
    ap = waveform(:,biggestChan);
    [~,mx] = max(ap);
    [~,mn] = min(ap);
    pyrWidths(i) = abs(mx-mn)/fs;
end

intWidths = nan(1,numel(intInds));
for i = 1:numel(intInds)
    waveform = neuronDataStruct(intInds(i)).waveforms;
    biggestChan = neuronDataStruct(intInds(i)).biggestChan;
    ap = waveform(:,biggestChan);
    [~,mx] = max(ap);
    [~,mn] = min(ap);
    intWidths(i) = abs(mx-mn)/fs;
end

[~,pyrMedIdx] = min(abs(pyrWidths - median(pyrWidths,'omitnan')));
[~,intMedIdx] = min(abs(intWidths - median(intWidths,'omitnan')));

pyrExampleInd = pyrInds(pyrMedIdx);
intExampleInd = intInds(intMedIdx);
end

function plotSingleWaveformInset(unitStruct, fs, waveColor, figTitle)
% plots one waveform with a 0.5 ms reference bar

waveforms = unitStruct.waveforms;
biggestChan = unitStruct.biggestChan;
ap = double(waveforms(:,biggestChan));

t_ms = (0:numel(ap)-1) / fs * 1000; % ms

figure('Color','w','Position',[100 100 260 220]);
ax = axes; hold(ax,'on');

plot(ax, t_ms, ap, 'Color', waveColor, 'LineWidth', 2);

% clean inset look
box(ax,'off');
ax.XColor = 'none';
ax.YColor = 'none';
ax.LineWidth = 1.2;

% add 0.5 ms scale bar
xRange = max(t_ms) - min(t_ms);
yRange = max(ap) - min(ap);

xBarStart = min(t_ms) + 0.08*xRange;
xBarEnd   = xBarStart + 0.5; % 0.5 ms
yBar      = min(ap) + 0.08*yRange;

plot(ax, [xBarStart xBarEnd], [yBar yBar], 'k', 'LineWidth', 2);
text(mean([xBarStart xBarEnd]), yBar - 0.08*yRange, '0.5 ms', ...
    'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize', 12);

title(ax, figTitle, 'FontSize', 14);

end
