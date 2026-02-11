function plotRasterWithBehaviorPatches(baseDir, t0, t1, labelType)
% plots raster snippet with behavior labels as background patches (no emg trace)
% inputs:
%   baseDir   : path to session ProcessedData folder
%   t0, t1    : snippet window in seconds
%   labelType : "umap" (only supported right now)

% example:
%   plotRasterWithBehaviorPatches("Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData", 600, 675, "umap")

arguments
    baseDir (1,1) string
    t0 (1,1) double
    t1 (1,1) double
    labelType (1,1) string = "umap"
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

%% ---- load behavior labels + sync mapping (only required vars!) ----
switch lower(labelType)
    case "umap"
        U = load(fullfile(baseDir,'UMAP.mat'), 'origDownsampEMGInd', 'regionAssignmentsFiltered');
        if ~isfield(U,'origDownsampEMGInd') || isempty(U.origDownsampEMGInd)
            error('UMAP.mat missing origDownsampEMGInd');
        end
        if ~isfield(U,'regionAssignmentsFiltered') || isempty(U.regionAssignmentsFiltered)
            error('UMAP.mat missing regionAssignmentsFiltered');
        end

        origDownsampEMGInd = double(U.origDownsampEMGInd(:)); % reduced->full emg idx
        regionRaw = double(U.regionAssignmentsFiltered(:)); % values 33..39
        regionCanon = regionRaw - 32; % map 33..39 -> 1..7

    otherwise
        error('labelType "%s" not supported yet. use "umap".', labelType);
end

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

neurInds_ms = origDownsampEMGInd * emgNeurSlope + emgNeurOffset;  % ms index in neural time base
labelTimes_s = neurInds_ms / 1000;                                  % seconds

% guard: sizes gotta agree
n = min(numel(labelTimes_s), numel(regionCanon));
labelTimes_s = labelTimes_s(1:n);
regionCanon = regionCanon(1:n);

%% ---- IMPORTANT: restrict labels to the plotting window BEFORE block finding ----
inWin = (labelTimes_s >= t0) & (labelTimes_s <= t1) & ~isnan(regionCanon);
if ~any(inWin)
    warning('no behavior labels found in this window. plotting raster only.');
end

tWin = labelTimes_s(inWin);
labWin = regionCanon(inWin);

% ensure row vectors for easier diff/block logic
tWin = tWin(:)';
labWin = labWin(:)';

%% ---- make figure ----
figure('Color','w','Position',[100 100 1000 500]);
ax = axes; hold on; %#ok<NASGU>

yl = [0, nNeur+1];

% 1) behavior patches (behind raster)
% stable colors 1..7
cmap = lines(7);

if ~isempty(labWin)
    % find contiguous blocks where label stays constant in the windowed vector
    changePts = [true, diff(labWin) ~= 0];
    blockStarts = find(changePts);
    blockStops  = [blockStarts(2:end)-1, numel(labWin)];

    for b = 1:numel(blockStarts)
        i0 = blockStarts(b);
        i1 = blockStops(b);

        beh = labWin(i0);
        if isnan(beh) || beh < 1 || beh > 7
            continue;
        end

        tStart = tWin(i0);
        tStop  = tWin(i1);

        % safety clip (should already be in window)
        tStart = max(tStart, t0);
        tStop  = min(tStop,  t1);

        if tStop > tStart
            patch([tStart tStart tStop tStop], [yl(1) yl(2) yl(2) yl(1)], cmap(beh,:), ...
                'EdgeColor','none', 'FaceAlpha',0.15);
        end
    end
end

% 2) raster on top
for iNeuron = 1:nNeur
    t = tsSec{iNeuron};
    t = t(t >= t0 & t <= t1);
    if ~isempty(t)
        plot(t, iNeuron*ones(size(t)), 'k.', 'MarkerSize', 2);
    end
end

xlim([t0 t1]);
ylim(yl);
xlabel('time (s)');
ylabel('neuron index');
title(sprintf('raster snippet with umap behavior patches | %.1f–%.1f s', t0, t1));
box off;

end
