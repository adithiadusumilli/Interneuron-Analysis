function plotExampleWaveformsForInset(baseDir)
% plots one example cortex pyramidal waveform and one example cortex interneuron waveform in separate small figures for use as insets

% also plots average cortex neuron-type percentages across animals:
%   - bars = mean % across animals
%   - dots = one dot per animal
%   - lines connect pyramidal and interneuron percentages within each animal

% j run: plotExampleWaveformsForInset("Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData")

arguments
    baseDir (1,1) string
end

%% ---- settings ----
conslidatedDataFoler = 'X:\David\AnalysesData';

animalFolders = {
    'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData', ...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData', ...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData', ...
    'X:\David\ArenaRecordings\D043-020525-ArenaRecording\ProcessedData'
};

fs = 30000; % hz

%% ---- load data for current animal ----
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

%% ---- pick one new example pyr and one new example int ----
pyrInds = cortexInds(cortexLabels == 0);
intInds = cortexInds(cortexLabels == 1);

if isempty(pyrInds) || isempty(intInds)
    error('could not find both pyramidal and interneuron cortex units.');
end

% choose more distinct examples:
% pyramidal = unit near upper quartile width
% interneuron = unit near lower quartile width
[pyrExampleInd, intExampleInd] = pickDistinctExamples(neuronDataStruct, pyrInds, intInds, fs);

%% ---- plot pyramidal waveform ----
plotSingleWaveformInset(neuronDataStruct(pyrExampleInd), fs, [0 0 1], 'Pyramidal Example');

%% ---- plot interneuron waveform ----
plotSingleWaveformInset(neuronDataStruct(intExampleInd), fs, [1 0 0], 'Interneuron Example');

%% ---- plot neuron-type percentages across animals ----
pyrPctAll = nan(numel(animalFolders),1);
intPctAll = nan(numel(animalFolders),1);

for a = 1:numel(animalFolders)
    thisBaseDir = string(animalFolders{a});

    Fr = load(fullfile(thisBaseDir,'NeuralFiringRates1msBins10msGauss.mat'),'cortexInds');
    thisCortexInds = double(Fr.cortexInds(:))';

    thisLabelsAll = double(classifications{a,1});
    thisCortexLabels = thisLabelsAll(thisCortexInds); % 0 = pyr, 1 = int

    nPyr = sum(thisCortexLabels == 0);
    nInt = sum(thisCortexLabels == 1);
    nTot = nPyr + nInt;

    if nTot > 0
        pyrPctAll(a) = 100 * nPyr / nTot;
        intPctAll(a) = 100 * nInt / nTot;
    end
end

meanPyrPct = mean(pyrPctAll, 'omitnan');
meanIntPct = mean(intPctAll, 'omitnan');

figure('Color','w','Position',[100 100 380 290]);
ax = axes; hold(ax,'on');

% average bars
b1 = bar(ax, 1, meanPyrPct, 'FaceColor', [0 0 1], 'EdgeColor', 'none');
b2 = bar(ax, 2, meanIntPct, 'FaceColor', [1 0 0], 'EdgeColor', 'none');
b1.FaceAlpha = 0.28;
b2.FaceAlpha = 0.28;

% connect each animal's pyr and int percentages
for a = 1:numel(animalFolders)
    plot(ax, [1 2], [pyrPctAll(a) intPctAll(a)], '-', ...
        'Color', [0.4 0.4 0.4], 'LineWidth', 1.1);
end

% dots for each animal
plot(ax, ones(size(pyrPctAll)), pyrPctAll, 'o', ...
    'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', 'k', 'MarkerSize', 6);

plot(ax, 2*ones(size(intPctAll)), intPctAll, 'o', ...
    'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'k', 'MarkerSize', 6);

xlim(ax, [0.4 2.6]);
xticks(ax, [1 2]);
xticklabels(ax, {'Pyramidal Neuron','Interneuron'});
ylabel(ax, 'Percent of M1 Units (%)');
title(ax, 'M1 Neuron-Type Percentages', 'FontSize', 16);

box(ax,'off');
ax.LineWidth = 1.2;
ax.FontSize = 16;
ax.TickDir = 'out';

end

function [pyrExampleInd, intExampleInd] = pickDistinctExamples(neuronDataStruct, pyrInds, intInds, fs)
% choose more distinct examples:
%   - pyramidal: closest to 75th percentile width
%   - interneuron: closest to 25th percentile width

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

pyrTarget = prctile(pyrWidths, 75);
intTarget = prctile(intWidths, 25);

[~,pyrIdx] = min(abs(pyrWidths - pyrTarget));
[~,intIdx] = min(abs(intWidths - intTarget));

pyrExampleInd = pyrInds(pyrIdx);
intExampleInd = intInds(intIdx);
end

function plotSingleWaveformInset(unitStruct, fs, waveColor, figTitle)
% plots one waveform emphasizing the peak-to-peak measurement with a 0.5 ms reference bar

waveforms = unitStruct.waveforms;
biggestChan = unitStruct.biggestChan;
ap = double(waveforms(:,biggestChan));

% make polarity consistent so trough comes before peak when possible
[~,mx0] = max(ap);
[~,mn0] = min(ap);
if mx0 < mn0
    ap = -ap;
end

% recompute after possible flip
[~,mx] = max(ap);
[~,mn] = min(ap);

t_ms = (0:numel(ap)-1) / fs * 1000; % ms

% crop around the peak-to-peak event
leftPad = 8;
rightPad = 12;
i0 = max(1, mn - leftPad);
i1 = min(numel(ap), mx + rightPad);

apCrop = ap(i0:i1);
tCrop_ms = t_ms(i0:i1);

% local indices for markers in cropped trace
mnLocal = mn - i0 + 1;
mxLocal = mx - i0 + 1;

figure('Color','w','Position',[100 100 260 220]);
ax = axes; hold(ax,'on');

plot(ax, tCrop_ms, apCrop, 'Color', waveColor, 'LineWidth', 2.5);

% mark trough and peak used for duration
plot(ax, tCrop_ms(mnLocal), apCrop(mnLocal), 'o', ...
    'MarkerFaceColor', waveColor, 'MarkerEdgeColor', waveColor, 'MarkerSize', 5);
plot(ax, tCrop_ms(mxLocal), apCrop(mxLocal), 'o', ...
    'MarkerFaceColor', waveColor, 'MarkerEdgeColor', waveColor, 'MarkerSize', 5);

% clean inset look
box(ax,'off');
ax.XColor = 'none';
ax.YColor = 'none';
ax.LineWidth = 1.2;

% tighten limits
xlim(ax, [min(tCrop_ms) max(tCrop_ms)]);
yPad = 0.18 * (max(apCrop) - min(apCrop));
ylim(ax, [min(apCrop)-yPad, max(apCrop)+0.10*yPad]);

% add 0.5 ms scale bar
xRange = max(tCrop_ms) - min(tCrop_ms);
yRange = diff(ylim(ax));

xBarStart = min(tCrop_ms) + 0.10*xRange;
xBarEnd   = xBarStart + 0.5; % 0.5 ms
yBar      = min(ylim(ax)) + 0.15*yRange;

plot(ax, [xBarStart xBarEnd], [yBar yBar], 'k', 'LineWidth', 2);
text(mean([xBarStart xBarEnd]), yBar - 0.10*yRange, '0.5 ms', ...
    'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize', 16);

title(ax, figTitle, 'FontSize', 16);

end
