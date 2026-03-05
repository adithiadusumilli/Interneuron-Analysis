function plotRasterWithBehaviorPatches(baseDir, t0, t1, labelType)
% plots a cortex-only raster snippet with behavior labels as background patches (no emg trace)

% inputs:
%   basedir: path to session ProcessedDdata folder
%   t0, t1: snippet window in seconds
%   labeltype: "classifier" (default), "manual", or "umap"
%
% run:
%   plotRasterWithBehaviorPatches("Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData", 600, 675, "classifier")

arguments
    baseDir (1,1) string
    t0 (1,1) double
    t1 (1,1) double
    labelType (1,1) string = "classifier"
end

fsNeur = 30000; % neuropixels hz

if t1 <= t0
    error('t1 must be > t0');
end

%% ---- load cortex indices from classifications ----
C = load(fullfile(baseDir,'AA_classifications.mat'));
cortexInds = local_getCortexInds(C);

if isempty(cortexInds)
    error('could not find cortex indices in AA_classifications.mat');
end

%% ---- load spikes ----
S = load(fullfile(baseDir,'neuronDataStruct.mat'),'neuronDataStruct');
neuronDataStruct = S.neuronDataStruct;

tsSec_all = cellfun(@(x) double(x)/fsNeur, {neuronDataStruct.timeStamps}, 'UniformOutput', false);
nAll = numel(tsSec_all);

% make cortexInds an index vector (not logical), and bounds-check
if islogical(cortexInds)
    cortexInds = find(cortexInds);
end
cortexInds = cortexInds(:)';
cortexInds = cortexInds(cortexInds >= 1 & cortexInds <= nAll);

if isempty(cortexInds)
    error('cortexInds was empty after bounds-checking against neuronDataStruct (%d units)', nAll);
end

% cortex-only spike times
tsSec = tsSec_all(cortexInds);
nNeur = numel(tsSec);

%% ---- load time-resolved behavior labels from umap + build canonical labels ----
U = load(fullfile(baseDir,'UMAP.mat'), ...
    'origDownsampEMGInd', ...
    'regionAssignmentsFiltered', ...
    'behvLabelsNoArt', ...
    'classifierLabels', ...
    'classifierBehvs', ...
    'analyzedBehaviors', ...
    'regionBehvAssignments');

% checks
req = {'origDownsampEMGInd','regionAssignmentsFiltered','behvLabelsNoArt','classifierLabels'};
for i = 1:numel(req)
    if ~isfield(U,req{i}) || isempty(U.(req{i}))
        error('UMAP.mat missing %s', req{i});
    end
end

origDownsampEMGInd = double(U.origDownsampEMGInd(:)); % reduced -> full emg idx
regionRaw = double(U.regionAssignmentsFiltered(:)); % raw region codes
manualRaw = double(U.behvLabelsNoArt(:)); % 0..nManual
classRaw = double(U.classifierLabels(:)); % 0..nClass

% canonical behavior name list and ordering (shared across animals)
manBehvNames = {'climbdown','climbup','eating','grooming', 'jumpdown','jumping','rearing','still','walkflat','walkgrid'};

% umap regions -> contiguous 1..7
regionCodes = unique(regionRaw(~isnan(regionRaw)));
nRegions = numel(regionCodes);
if nRegions ~= 7
    error('expected 7 umap regions, found %d', nRegions);
end

umapCanon = nan(size(regionRaw)); % 1..7
for r = 1:nRegions
    umapCanon(regionRaw == regionCodes(r)) = r;
end

% manual labels -> canonical 0..10
manualCanon = zeros(size(manualRaw)); % default 0 = unlabeled
if isfield(U,'analyzedBehaviors') && ~isempty(U.analyzedBehaviors)
    analyzedBehaviors = U.analyzedBehaviors;

    manBehvNumbers = zeros(1, numel(analyzedBehaviors)); % per-animal -> canonical (0..10)
    for iBehv = 1:numel(analyzedBehaviors)
        idx = find(strcmp(analyzedBehaviors{iBehv}, manBehvNames), 1);
        manBehvNumbers(iBehv) = iff(isempty(idx), 0, idx);
    end

    for i = 1:numel(manualRaw)
        v = manualRaw(i);
        if v == 0
            manualCanon(i) = 0;
        else
            if v >= 1 && v <= numel(manBehvNumbers)
                manualCanon(i) = manBehvNumbers(v);
            else
                manualCanon(i) = 0;
            end
        end
    end
end

% classifier labels -> canonical 0..10
classifierCanon = zeros(size(classRaw)); % default 0 = unlabeled
if isfield(U,'classifierBehvs') && ~isempty(U.classifierBehvs)
    classifierBehvs = U.classifierBehvs;

    classBehvNumbers = zeros(1, numel(classifierBehvs)); % per-animal -> canonical (0..10)
    for iBehv = 1:numel(classifierBehvs)
        idx = find(strcmp(classifierBehvs{iBehv}, manBehvNames), 1);
        classBehvNumbers(iBehv) = iff(isempty(idx), 0, idx);
    end

    for i = 1:numel(classRaw)
        v = classRaw(i);
        if v == 0
            classifierCanon(i) = 0;
        else
            if v >= 1 && v <= numel(classBehvNumbers)
                classifierCanon(i) = classBehvNumbers(v);
            else
                classifierCanon(i) = 0;
            end
        end
    end
end

% choose label stream to plot
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
               (round(frameEMGSamples{1}{end}(end)/20) - round(frameEMGSamples{1}{1}(1)/20));
emgNeurOffset = round(frameNeuropixelSamples{1}{1}(1)/30) - emgNeurSlope*round(frameEMGSamples{1}{1}(1)/20);

neurInds_ms = origDownsampEMGInd * emgNeurSlope + emgNeurOffset; % ms index in neural time base
labelTimes_s = neurInds_ms / 1000; % seconds

% guard: sizes gotta agree across time and labels
n = min(numel(labelTimes_s), numel(labelsCanon));
labelTimes_s = labelTimes_s(1:n);
labelsCanon = labelsCanon(1:n);

%% ---- restrict labels to plotting window before block finding ----
inWin = (labelTimes_s >= t0) & (labelTimes_s <= t1) & ~isnan(labelsCanon);
tWin = labelTimes_s(inWin);
labWin = labelsCanon(inWin);

tWin = tWin(:)';   % row vectors for diff logic
labWin = labWin(:)';

%% ---- make figure ----
fig = figure('Color','w','Position',[100 100 1100 520]);
ax = axes(fig);
hold(ax,'on');

% behavior patches behind raster
yl = [0, nNeur+1];
cmap = local_behavior_cmap(nColors);

if ~isempty(labWin)
    changePts = [true, diff(labWin) ~= 0];
    blockStarts = find(changePts);
    blockStops = [blockStarts(2:end)-1, numel(labWin)];

    for b = 1:numel(blockStarts)
        i0 = blockStarts(b);
        i1 = blockStops(b);

        beh = labWin(i0);
        cInd = local_label_to_color_index(beh, nColors);
        if isnan(cInd)
            continue;
        end

        tStart = max(tWin(i0), t0);
        tStop = min(tWin(i1), t1);

        if tStop > tStart
            patch(ax, [tStart tStart tStop tStop], [yl(1) yl(2) yl(2) yl(1)], cmap(cInd,:), ...
                'EdgeColor','none', 'FaceAlpha',0.12);
        end
    end
end

% raster on top (publication-style: vertical tick marks, vectorized per neuron)
tickHalfHeight = 0.35;

for iNeuron = 1:nNeur
    t = tsSec{iNeuron};
    t = t(t >= t0 & t <= t1);
    if isempty(t)
        continue;
    end

    x = [t(:)'; t(:)'];
    y = [(iNeuron - tickHalfHeight)*ones(1,numel(t)); (iNeuron + tickHalfHeight)*ones(1,numel(t))];
    line(ax, x, y, 'Color', 'k', 'LineWidth', 0.5);
end

% axes formatting
xlim(ax, [t0 t1]);
ylim(ax, yl);

ax.FontSize  = 14;
ax.TickDir   = 'out';
ax.LineWidth = 1;

xlabel(ax, 'time (s)', 'FontSize', 16);
ylabel(ax, 'cortex neuron index', 'FontSize', 16);
title(ax, sprintf('cortex raster with %s behavior patches | %.1f–%.1f s | n=%d units', ...
    labelType, t0, t1, nNeur), 'FontSize', 18);

box(ax, 'off');

% legend
behNames = local_behavior_names(U, manBehvNames, nColors);
h = gobjects(1, numel(behNames));
for k = 1:numel(behNames)
    h(k) = patch(ax, nan, nan, cmap(k,:), 'EdgeColor','none', 'FaceAlpha',0.12);
end
legend(ax, h, cellstr(behNames), 'Location', 'eastoutside');

end

%% ---------------- local helpers ----------------
function cortexInds = local_getCortexInds(C)
% tries common field names for cortex indices

cortexInds = [];

candidates = {'cortexInds','cortex_inds','cortexUnits','cortexUnitInds','cortexIdx','cortexInd'};
for i = 1:numel(candidates)
    f = candidates{i};
    if isfield(C,f) && ~isempty(C.(f))
        cortexInds = C.(f);
        return;
    end
end

% sometimes saved as a struct inside the mat
if isfield(C,'classifications') && isstruct(C.classifications)
    if isfield(C.classifications,'cortexInds') && ~isempty(C.classifications.cortexInds)
        cortexInds = C.classifications.cortexInds;
        return;
    end
end
end

function cmap = local_behavior_cmap(nColors)
% high-contrast fixed palette
% for nColors==11: 1 is "unlabeled" (light gray), 2..11 are canonical behaviors 1..10
% for nColors==7:  7 distinct colors

if nColors == 11
    cmap = [
        0.85 0.85 0.85;  % unlabeled
        0.00 0.45 0.70;  % climbdown
        0.90 0.62 0.00;  % climbup
        0.00 0.62 0.45;  % eating
        0.80 0.47 0.65;  % grooming
        0.34 0.71 0.91;  % jumpdown
        0.84 0.37 0.00;  % jumping
        0.94 0.89 0.26;  % rearing
        0.49 0.18 0.56;  % still
        0.20 0.29 0.37;  % walkflat
        0.64 0.08 0.18;  % walkgrid
    ];
else
    cmap = [
        0.00 0.45 0.70;
        0.90 0.62 0.00;
        0.00 0.62 0.45;
        0.80 0.47 0.65;
        0.34 0.71 0.91;
        0.84 0.37 0.00;
        0.49 0.18 0.56;
    ];
end
end

function cInd = local_label_to_color_index(beh, nColors)
% converts canonical label values to colormap row indices
% returns nan if label is invalid

cInd = nan;

if nColors == 7
    if isnan(beh) || beh < 1 || beh > 7
        return;
    end
    cInd = beh;
else
    if isnan(beh) || beh < 0 || beh > 10
        return;
    end
    cInd = beh + 1; % 0..10 -> 1..11
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
if cond
    out = a;
else
    out = b;
end
end
