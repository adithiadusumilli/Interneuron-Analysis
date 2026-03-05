function plotRasterAndTwoEMGWithBehaviorPatches(baseDir, t0, t1, emgChan1, emgChan2, labelType)
% make 3-panel vertical plot:
%   1) cortex-only raster (top)
%   2) emg channel emgChan1 (middle)
%   3) emg channel emgChan2 (bottom)
% all panels share the same behavior background patches for the same time window

% j run: plotRasterAndTwoEMGWithBehaviorPatches("Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData", 600, 675, 1, 2, "classifier")

arguments
    baseDir (1,1) string
    t0 (1,1) double
    t1 (1,1) double
    emgChan1 (1,1) double = 1
    emgChan2 (1,1) double = 2
    labelType (1,1) string = "classifier"
end

if t1 <= t0
    error('t1 must be > t0');
end

fsNeur = 30000; % neuropixels hz
fsEmg  = 1000;  % downsamp emg hz

%% ---- load cortex indices from firing rate file ----
F = load(fullfile(baseDir,'NeuralFiringRates1msBins10msGauss.mat'),'cortexInds');
if ~isfield(F,'cortexInds') || isempty(F.cortexInds)
    error('could not find cortexInds in NeuralFiringRates1msBins10msGauss.mat');
end
cortexInds = F.cortexInds;

%% ---- load spikes ----
S = load(fullfile(baseDir,'neuronDataStruct.mat'),'neuronDataStruct');
neuronDataStruct = S.neuronDataStruct;

tsSec_all = cellfun(@(x) double(x)/fsNeur, {neuronDataStruct.timeStamps}, 'UniformOutput', false);
nAll = numel(tsSec_all);

if islogical(cortexInds)
    cortexInds = find(cortexInds);
end
cortexInds = cortexInds(:)';
cortexInds = cortexInds(cortexInds >= 1 & cortexInds <= nAll);

if isempty(cortexInds)
    error('cortexInds was empty after bounds-checking against neuronDataStruct (%d units)', nAll);
end

tsSec = tsSec_all(cortexInds);
nNeur = numel(tsSec);

%% ---- load emg once (both channels will index into this) ----
E = load(fullfile(baseDir,'EMG1ms.mat'));
if isfield(E,'downsampEMG')
    emgAll = E.downsampEMG;
elseif isfield(E,'EMG')
    emgAll = E.EMG;
else
    error('could not find emg variable in EMG1ms.mat');
end
emgTime = (0:size(emgAll,2)-1) / fsEmg; % assumes channels x time if matrix

%% ---- load time-resolved behavior labels from umap + build canonical labels ----
U = load(fullfile(baseDir,'UMAP.mat'), ...
    'origDownsampEMGInd', ...
    'regionAssignmentsFiltered', ...
    'behvLabelsNoArt', ...
    'classifierLabels', ...
    'classifierBehvs', ...
    'analyzedBehaviors', ...
    'regionBehvAssignments');

req = {'origDownsampEMGInd','regionAssignmentsFiltered','behvLabelsNoArt','classifierLabels'};
for i = 1:numel(req)
    if ~isfield(U,req{i}) || isempty(U.(req{i}))
        error('UMAP.mat missing %s', req{i});
    end
end

origDownsampEMGInd = double(U.origDownsampEMGInd(:));
regionRaw = double(U.regionAssignmentsFiltered(:));
manualRaw = double(U.behvLabelsNoArt(:));
classRaw  = double(U.classifierLabels(:));

manBehvNames = {'climbdown','climbup','eating','grooming', 'jumpdown','jumping','rearing','still','walkflat','walkgrid'};

% umap regions -> contiguous 1..7
regionCodes = unique(regionRaw(~isnan(regionRaw)));
nRegions = numel(regionCodes);
if nRegions ~= 7
    error('expected 7 umap regions, found %d', nRegions);
end
umapCanon = nan(size(regionRaw));
for r = 1:nRegions
    umapCanon(regionRaw == regionCodes(r)) = r;
end

% manual labels -> canonical 0..10
manualCanon = zeros(size(manualRaw));
if isfield(U,'analyzedBehaviors') && ~isempty(U.analyzedBehaviors)
    analyzedBehaviors = U.analyzedBehaviors;
    manBehvNumbers = zeros(1, numel(analyzedBehaviors));
    for iBehv = 1:numel(analyzedBehaviors)
        idx = find(strcmp(analyzedBehaviors{iBehv}, manBehvNames), 1);
        manBehvNumbers(iBehv) = iff(isempty(idx), 0, idx);
    end
    for i0 = 1:numel(manualRaw)
        v = manualRaw(i0);
        if v == 0
            manualCanon(i0) = 0;
        elseif v >= 1 && v <= numel(manBehvNumbers)
            manualCanon(i0) = manBehvNumbers(v);
        else
            manualCanon(i0) = 0;
        end
    end
end

% classifier labels -> canonical 0..10
classifierCanon = zeros(size(classRaw));
if isfield(U,'classifierBehvs') && ~isempty(U.classifierBehvs)
    classifierBehvs = U.classifierBehvs;
    classBehvNumbers = zeros(1, numel(classifierBehvs));
    for iBehv = 1:numel(classifierBehvs)
        idx = find(strcmp(classifierBehvs{iBehv}, manBehvNames), 1);
        classBehvNumbers(iBehv) = iff(isempty(idx), 0, idx);
    end
    for i0 = 1:numel(classRaw)
        v = classRaw(i0);
        if v == 0
            classifierCanon(i0) = 0;
        elseif v >= 1 && v <= numel(classBehvNumbers)
            classifierCanon(i0) = classBehvNumbers(v);
        else
            classifierCanon(i0) = 0;
        end
    end
end

% choose label stream
switch lower(labelType)
    case "umap"
        labelsCanon = umapCanon; % 1..7
        nColors = 7;
    case "manual"
        labelsCanon = manualCanon; % 0..10
        nColors = 11;
    case "classifier"
        labelsCanon = classifierCanon; % 0..10
        nColors = 11;
    otherwise
        error('labeltype "%s" not supported. use "classifier", "manual", or "umap".', labelType);
end

%% ---- load sync mapping ----
V = load(fullfile(baseDir,'VideoSyncFrames.mat'), 'frameNeuropixelSamples', 'frameEMGSamples');
if ~isfield(V,'frameNeuropixelSamples') || ~isfield(V,'frameEMGSamples')
    error('VideoSyncFrames.mat missing frameNeuropixelSamples and/or frameEMGSamples');
end
frameNeuropixelSamples = V.frameNeuropixelSamples;
frameEMGSamples = V.frameEMGSamples;

%% ---- map reduced umap index -> neuropixels time (seconds) ----
emgNeurSlope = (round(frameNeuropixelSamples{1}{end}(end)/30) - round(frameNeuropixelSamples{1}{1}(1)/30)) / ...
               (round(frameEMGSamples{1}{end}(end)/20)     - round(frameEMGSamples{1}{1}(1)/20));
emgNeurOffset = round(frameNeuropixelSamples{1}{1}(1)/30) - emgNeurSlope*round(frameEMGSamples{1}{1}(1)/20);

neurInds_ms  = origDownsampEMGInd * emgNeurSlope + emgNeurOffset;
labelTimes_s = neurInds_ms / 1000;

n = min(numel(labelTimes_s), numel(labelsCanon));
labelTimes_s = labelTimes_s(1:n);
labelsCanon  = labelsCanon(1:n);

% restrict to plotting window + remove nan labels
inWin = (labelTimes_s >= t0) & (labelTimes_s <= t1) & ~isnan(labelsCanon);
tWin  = labelTimes_s(inWin);
labWin = labelsCanon(inWin);
tWin = tWin(:)'; labWin = labWin(:)';

% precompute behavior blocks once (tStart, tStop, colorIndex)
blocks = local_make_blocks(tWin, labWin, t0, t1, nColors);

% shared colormap + names for legend
cmap = local_behavior_cmap(nColors);
behNames = local_behavior_names(U, manBehvNames, nColors);

%% ---- make figure + axes ----
fig = figure('Color','w','Position',[100 100 1100 900]);
tl = tiledlayout(fig, 3, 1, 'TileSpacing','compact', 'Padding','compact');

ax1 = nexttile(tl, 1); hold(ax1,'on'); % raster
ax2 = nexttile(tl, 2); hold(ax2,'on'); % emg 1
ax3 = nexttile(tl, 3); hold(ax3,'on'); % emg 2

%% ---- panel 1: raster ----
ylRaster = [0, nNeur+1];
local_draw_patches(ax1, blocks, ylRaster, cmap);

tickHalfHeight = 0.35;
for iNeuron = 1:nNeur
    t = tsSec{iNeuron};
    t = t(t >= t0 & t <= t1);
    if isempty(t), continue; end
    x = [t(:)'; t(:)'];
    y = [(iNeuron - tickHalfHeight)*ones(1,numel(t)); (iNeuron + tickHalfHeight)*ones(1,numel(t))];
    line(ax1, x, y, 'Color', 'k', 'LineWidth', 0.5);
end

xlim(ax1, [t0 t1]);
ylim(ax1, ylRaster);
ax1.FontSize = 14; ax1.TickDir = 'out'; ax1.LineWidth = 1;
ylabel(ax1, 'M1 Neuron Index', 'FontSize', 16);
title(ax1, sprintf('M1 Raster + EMG | %s labels | %.1f–%.1f s | n=%d units', ...
    labelType, t0, t1, nNeur), 'FontSize', 18);
box(ax1,'off');

%% ---- panel 2: emg channel 1 ----
local_plot_emg_panel(ax2, emgAll, emgTime, emgChan1, t0, t1, blocks, cmap);
ylabel(ax2, sprintf('EMG ch %d (a.u.)', emgChan1), 'FontSize', 16);
box(ax2,'off');

%% ---- panel 3: emg channel 2 ----
local_plot_emg_panel(ax3, emgAll, emgTime, emgChan2, t0, t1, blocks, cmap);
ylabel(ax3, sprintf('EMG ch %d (a.u.)', emgChan2), 'FontSize', 16);
xlabel(ax3, 'Time (s)', 'FontSize', 16);
box(ax3,'off');

%% ---- align x across panels ----
linkaxes([ax1 ax2 ax3], 'x');

%% ---- one shared legend for behavior colors (right side) ----
h = gobjects(1, numel(behNames));
for k = 1:numel(behNames)
    h(k) = patch(ax1, nan, nan, cmap(k,:), 'EdgeColor','none', 'FaceAlpha',0.30);
end
legend(ax1, h, cellstr(behNames), 'Location', 'eastoutside');

end

%% ---------------- local helpers ----------------
function blocks = local_make_blocks(tWin, labWin, t0, t1, nColors)
% returns struct array with fields: tStart, tStop, cInd
blocks = struct('tStart',{},'tStop',{},'cInd',{});

if isempty(labWin), return; end

changePts = [true, diff(labWin) ~= 0];
blockStarts = find(changePts);
blockStops  = [blockStarts(2:end)-1, numel(labWin)];

bCount = 0;
for b = 1:numel(blockStarts)
    i0 = blockStarts(b);
    i1 = blockStops(b);

    beh = labWin(i0);
    cInd = local_label_to_color_index(beh, nColors);
    if isnan(cInd), continue; end

    tStart = max(tWin(i0), t0);
    tStop  = min(tWin(i1), t1);
    if tStop <= tStart, continue; end

    bCount = bCount + 1;
    blocks(bCount).tStart = tStart;
    blocks(bCount).tStop  = tStop;
    blocks(bCount).cInd   = cInd;
end
end

function local_draw_patches(ax, blocks, yl, cmap)
% draws all behavior patches for a given axis/y-limits
for b = 1:numel(blocks)
    c = cmap(blocks(b).cInd,:);
    patch(ax, [blocks(b).tStart blocks(b).tStart blocks(b).tStop blocks(b).tStop], ...
              [yl(1) yl(2) yl(2) yl(1)], c, 'EdgeColor','none', 'FaceAlpha',0.30);
end
end

function local_plot_emg_panel(ax, emgAll, emgTime, emgChannel, t0, t1, blocks, cmap)
% plots selected emg channel + draws shared patches behind it
if isvector(emgAll)
    if emgChannel ~= 1
        error('emg is a vector in this file; only emgChannel=1 is valid.');
    end
    emg = emgAll(:);
else
    if emgChannel < 1 || emgChannel > size(emgAll,1)
        error('emgChannel %d out of range (1..%d)', emgChannel, size(emgAll,1));
    end
    emg = emgAll(emgChannel, :)';
end

keep = (emgTime >= t0) & (emgTime <= t1);
hEmg = plot(ax, emgTime(keep), emg(keep), 'k', 'LineWidth', 1);

xlim(ax, [t0 t1]);
yl = get(ax,'YLim');

local_draw_patches(ax, blocks, yl, cmap);
uistack(hEmg,'top');

ax.FontSize = 14;
ax.TickDir = 'out';
ax.LineWidth = 1;
end

function cmap = local_behavior_cmap(nColors)
if nColors == 11
    cmap = [
        0.85 0.85 0.85;  % unlabeled
        0.00 0.45 0.85;  % climbdown
        0.60 0.00 0.00;  % climbup
        0.00 0.62 0.45;  % eating
        0.93 0.69 0.13;  % grooming
        0.49 0.18 0.56;  % jumpdown
        0.30 0.75 0.93;  % jumping
        0.20 0.90 0.20;  % rearing
        0.90 0.00 0.60;  % still
        0.00 0.00 0.00;  % walkflat
        0.80 0.35 0.00   % walkgrid
    ];
else
    cmap = [
        0.00 0.45 0.74
        0.85 0.33 0.10
        0.00 0.62 0.45
        0.93 0.69 0.13
        0.49 0.18 0.56
        0.30 0.75 0.93
        1.00 0.50 0.00
    ];
end
end

function cInd = local_label_to_color_index(beh, nColors)
cInd = nan;
if nColors == 7
    if isnan(beh) || beh < 1 || beh > 7, return; end
    cInd = beh;
else
    if isnan(beh) || beh < 0 || beh > 10, return; end
    cInd = beh + 1;
end
end

function behNames = local_behavior_names(U, manBehvNames, nColors)
if nColors == 7
    if isfield(U,'regionBehvAssignments') && ~isempty(U.regionBehvAssignments)
        tmp = string(U.regionBehvAssignments(:)');
        tmp = tmp(1:min(7,numel(tmp)));
        if numel(tmp) < 7
            tmp(end+1:7) = "umap region " + (numel(tmp)+1:7);
        end
        behNames = tmp;
    else
        behNames = "umap region " + (1:7);
    end
else
    behNames = ["unlabeled", string(manBehvNames)];
end
end

function out = iff(cond, a, b)
if cond, out = a; else, out = b; end
end
