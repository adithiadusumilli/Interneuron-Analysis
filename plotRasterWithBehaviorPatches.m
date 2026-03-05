function plotRasterWithBehaviorPatches(baseDir, t0, t1, labelType)
% plots raster snippet with behavior labels as background patches (no emg trace)
% inputs:
%   basedir: path to session processedddata folder
%   t0, t1: snippet window in seconds
%   labeltype: "classifier" (default), "manual", or "umap"

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

%% ---- load spikes ----
S = load(fullfile(baseDir,'neuronDataStruct.mat'),'neuronDataStruct');
neuronDataStruct = S.neuronDataStruct;

tsSec = cellfun(@(x) double(x)/fsNeur, {neuronDataStruct.timeStamps}, 'UniformOutput', false);
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
% indices 1..10 in canonical space always correspond to these names
manBehvNames = {'climbdown','climbup','eating','grooming', ...
                'jumpdown','jumping','rearing','still','walkflat','walkgrid'};

%% ---- umap regions -> contiguous 1..7 ----
regionCodes = unique(regionRaw(~isnan(regionRaw)));
nRegions = numel(regionCodes);
if nRegions ~= 7
    error('expected 7 umap regions, found %d', nRegions);
end

umapCanon = nan(size(regionRaw)); % 1..7
for r = 1:nRegions
    umapCanon(regionRaw == regionCodes(r)) = r;
end

%% ---- manual labels -> canonical 0..10 ----
manualCanon = zeros(size(manualRaw)); % default 0 = unlabeled
if isfield(U,'analyzedBehaviors') && ~isempty(U.analyzedBehaviors)
    analyzedBehaviors = U.analyzedBehaviors;

    manBehvNumbers = zeros(1, numel(analyzedBehaviors)); % per-animal -> canonical (0..10)
    for iBehv = 1:numel(analyzedBehaviors)
        idx = find(strcmp(analyzedBehaviors{iBehv}, manBehvNames), 1);
        if isempty(idx)
            manBehvNumbers(iBehv) = 0;
        else
            manBehvNumbers(iBehv) = idx;
        end
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
else
    warning('UMAP.mat missing analyzedBehaviors; manual canonical labels will be 0/unlabeled.');
end

%% ---- classifier labels -> canonical 0..10 ----
classifierCanon = zeros(size(classRaw)); % default 0 = unlabeled
if isfield(U,'classifierBehvs') && ~isempty(U.classifierBehvs)
    classifierBehvs = U.classifierBehvs;

    classBehvNumbers = zeros(1, numel(classifierBehvs)); % per-animal -> canonical (0..10)
    for iBehv = 1:numel(classifierBehvs)
        idx = find(strcmp(classifierBehvs{iBehv}, manBehvNames), 1);
        if isempty(idx)
            classBehvNumbers(iBehv) = 0;
        else
            classBehvNumbers(iBehv) = idx;
        end
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
else
    warning('UMAP.mat missing classifierBehvs; classifier canonical labels will be 0/unlabeled.');
end

%% ---- choose label stream to plot ----
switch lower(labelType)
    case "umap"
        labelsCanon = umapCanon;       % 1..7
        nColors = 7;
    case "manual"
        labelsCanon = manualCanon;     % 0..10
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
% this matches original mapping logic
emgNeurSlope = (round(frameNeuropixelSamples{1}{end}(end)/30) - round(frameNeuropixelSamples{1}{1}(1)/30)) / ...
               (round(frameEMGSamples{1}{end}(end)/20)     - round(frameEMGSamples{1}{1}(1)/20));
emgNeurOffset = round(frameNeuropixelSamples{1}{1}(1)/30) - emgNeurSlope*round(frameEMGSamples{1}{1}(1)/20);

neurInds_ms = origDownsampEMGInd * emgNeurSlope + emgNeurOffset; % ms index in neural time base
labelTimes_s = neurInds_ms / 1000; % seconds

% guard: sizes gotta agree across time and labels
n = min(numel(labelTimes_s), numel(labelsCanon));
labelTimes_s = labelTimes_s(1:n);
labelsCanon  = labelsCanon(1:n);

%% ---- restrict labels to plotting window before block finding ----
inWin = (labelTimes_s >= t0) & (labelTimes_s <= t1) & ~isnan(labelsCanon);
if ~any(inWin)
    warning('no behavior labels found in this window. plotting raster only.');
end

tWin   = labelTimes_s(inWin);
labWin = labelsCanon(inWin);

% ensure row vectors for easier diff/block logic
tWin   = tWin(:)';
labWin = labWin(:)';

%% ---- make figure ----
fig = figure('Color','w','Position',[100 100 1000 500]);
ax = axes(fig);
hold(ax,'on');

% --- publication-style tweak: use vectorized raster with line objects (much cleaner than many plot calls)
% this makes the raster sharper + faster, and you can control spike "tick" height
tickHalfHeight = 0.35;  % controls how tall each spike tick is

% behavior patches behind raster
yl = [0, nNeur+1];
cmap = lines(nColors);

if ~isempty(labWin)
    changePts = [true, diff(labWin) ~= 0];
    blockStarts = find(changePts);
    blockStops  = [blockStarts(2:end)-1, numel(labWin)];

    for b = 1:numel(blockStarts)
        i0 = blockStarts(b);
        i1 = blockStops(b);

        beh = labWin(i0);

        if nColors == 7
            if isnan(beh) || beh < 1 || beh > 7
                continue;
            end
            cInd = beh;
        else
            if isnan(beh) || beh < 0 || beh > 10
                continue;
            end
            cInd = beh + 1; % 0..10 -> 1..11
        end

        tStart = max(tWin(i0), t0);
        tStop  = min(tWin(i1), t1);

        if tStop > tStart
            patch(ax, [tStart tStart tStop tStop], [yl(1) yl(2) yl(2) yl(1)], cmap(cInd,:), ...
                'EdgeColor','none', 'FaceAlpha',0.12);
        end
    end
end

% raster on top: vectorize each neuron as vertical ticks (fast + clean)
for iNeuron = 1:nNeur
    t = tsSec{iNeuron};
    t = t(t >= t0 & t <= t1);
    if isempty(t)
        continue;
    end

    % build a 2-by-n matrix for x and y so line() draws many vertical ticks at once
    x = [t(:)'; t(:)'];
    y = [(iNeuron - tickHalfHeight)*ones(1,numel(t)); (iNeuron + tickHalfHeight)*ones(1,numel(t))];

    line(ax, x, y, 'Color', 'k', 'LineWidth', 0.5);
end

% axes formatting
xlim(ax, [t0 t1]);
ylim(ax, yl);

ax.FontSize = 14; % tick label size
ax.TickDir = 'out'; % cleaner
ax.LineWidth = 1; % thicker axes
% ax.YDir = 'reverse'; 

xlabel(ax, 'time (s)', 'FontSize', 16);
ylabel(ax, 'neuron index', 'FontSize', 16);
title(ax, sprintf('raster snippet with %s behavior patches | %.1f–%.1f s', labelType, t0, t1), ...
    'FontSize', 18);

box(ax, 'off');

%% ---- legend mapping ----
if nColors == 7
    if isfield(U,'regionBehvAssignments') && ~isempty(U.regionBehvAssignments)
        behNames = string(U.regionBehvAssignments(:)');
        behNames = behNames(1:min(7,numel(behNames)));
        if numel(behNames) < 7
            behNames(end+1:7) = "umap region " + (numel(behNames)+1:7);
        end
    else
        behNames = "umap region " + (1:7);
    end
else
    behNames = ["unlabeled", string(manBehvNames)];
end

h = gobjects(1, numel(behNames));
for k = 1:numel(behNames)
    h(k) = patch(ax, nan, nan, cmap(k,:), 'EdgeColor','none', 'FaceAlpha',0.12);
end
legend(ax, h, cellstr(behNames), 'Location', 'eastoutside');

end
