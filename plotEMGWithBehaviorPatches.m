function plotEMGWithBehaviorPatches(baseDir, t0, t1, emgChannel, labelType)
% plotemgwithbehaviorpatches

% plots an emg snippet with behavior labels as background patches

% inputs:
%   basedir : session processedddata folder
%   t0, t1 : window in seconds
%   emgchannel : which channel to plot (since multichannel) 
%   labeltype : "classifier" (default), "manual", or "umap"

% j run plotEMGWithBehaviorPatches("Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData", 600, 675, 1, "classifier")

arguments
    baseDir (1,1) string
    t0 (1,1) double
    t1 (1,1) double
    emgChannel (1,1) double = 1
    labelType (1,1) string = "classifier"
end

fsEmg = 1000; % bc downsamp from 20khz

if t1 <= t0
    error('t1 must be > t0');
end

%% ---- load emg ----
E = load(fullfile(baseDir,'EMG1ms.mat'));
if isfield(E,'downsampEMG')
    emgAll = E.downsampEMG;
elseif isfield(E,'EMG')
    emgAll = E.EMG;
else
    error('could not find emg variable in EMG1ms.mat');
end

% handle channel selection
if isvector(emgAll)
    emg = emgAll(:);
else
    if emgChannel < 1 || emgChannel > size(emgAll,2)
        error('emgChannel %d out of range (1..%d)', emgChannel, size(emgAll,2));
    end
    emg = emgAll(:, emgChannel);
end

emgTime = (0:numel(emg)-1) / fsEmg;

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

%% ---- map reduced umap index -> emg time (seconds) ----
% note: we only need label times in seconds; same slope/offset logic as before
emgNeurSlope = (round(frameNeuropixelSamples{1}{end}(end)/30) - round(frameNeuropixelSamples{1}{1}(1)/30)) / ...
               (round(frameEMGSamples{1}{end}(end)/20) - round(frameEMGSamples{1}{1}(1)/20));
emgNeurOffset = round(frameNeuropixelSamples{1}{1}(1)/30) - emgNeurSlope*round(frameEMGSamples{1}{1}(1)/20);

neurInds_ms = origDownsampEMGInd * emgNeurSlope + emgNeurOffset; % ms in neural time base
labelTimes_s = double(neurInds_ms) / 1000; % seconds

% guard: sizes gotta agree
n = min(numel(labelTimes_s), numel(labelsCanon));
labelTimes_s = labelTimes_s(1:n);
labelsCanon = labelsCanon(1:n);

%% ---- restrict labels to plotting window before block finding ----
inWin = (labelTimes_s >= t0) & (labelTimes_s <= t1) & ~isnan(labelsCanon);
tWin = labelTimes_s(inWin);
labWin = labelsCanon(inWin);

% ensure row vectors for diff/block logic
tWin = tWin(:)';
labWin = labWin(:)';

%% ---- make figure ----
fig = figure('Color','w','Position',[100 100 1000 350]);
ax = axes(fig);
hold(ax,'on');

% plot emg first so we have y-limits
keep = (emgTime >= t0) & (emgTime <= t1);
plot(ax, emgTime(keep), emg(keep), 'k', 'LineWidth', 1);

xlim(ax, [t0 t1]);
yl = get(ax, 'YLim');

% behavior patches behind emg (same contiguous-block logic as raster)
% use a fixed high-contrast palette (better separation than lines())
if nColors == 11
    cmap = [
        0.85 0.85 0.85;  % unlabeled (light gray)

        0.00 0.45 0.85;  % climbdown (blue)
        0.60 0.00 0.00;  % climbup
        0.00 0.62 0.45;  % eating (green)
        0.93 0.69 0.13;  % grooming (yellow)

        0.49 0.18 0.56;  % jumpdown (purple)
        0.30 0.75 0.93;  % jumping (cyan)
        0.20 0.90 0.20;  % rearing (lime green)

        0.90 0.00 0.60;  % still (magenta)
        0.00 0.00 0.00;  % walkflat (black)

        0.80 0.35 0.00   % walkgrid (dark orange-red)
    ];
else
    % UMAP case
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

if ~isempty(labWin)
    changePts = [true, diff(labWin) ~= 0];
    blockStarts = find(changePts);
    blockStops = [blockStarts(2:end)-1, numel(labWin)];

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
        tStop = min(tWin(i1), t1);

        if tStop > tStart
            patch(ax, [tStart tStart tStop tStop], [yl(1) yl(2) yl(2) yl(1)], cmap(cInd,:), ...
                'EdgeColor','none', 'FaceAlpha',0.12);
        end
    end
end

% bring emg to front
uistack(findobj(ax,'Type','line'),'top');

% axes formatting (same vibe as raster)
ax.FontSize = 14;
ax.TickDir = 'out';
ax.LineWidth = 1;

xlabel(ax, 'time (s)', 'FontSize', 16);
ylabel(ax, 'emg (a.u.)', 'FontSize', 16);
title(ax, sprintf('emg snippet with %s behavior patches | %.1f–%.1f s', labelType, t0, t1), ...
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
