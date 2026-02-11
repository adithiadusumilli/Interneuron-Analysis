function plotRasterWithBehaviorPatches(baseDir, t0, t1, labelType)
% plotRasterWithBehaviorPatches

% plots raster snippet with behavior labels as background patches
% (no emg trace shown)

% inputs:
%   baseDir   : path to session ProcessedData folder
%   t0, t1    : snippet window in seconds
%   labelType : 'umap' (default). you can extend later for manual labels.

% example:
%   plotRasterWithBehaviorPatches('Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData',600,675,'umap')

arguments
    baseDir (1,1) string
    t0 (1,1) double
    t1 (1,1) double
    labelType (1,1) string = "umap"
end

fsNeur = 30000;  % neuropixels hz

%% ---- load spikes ----
S = load(fullfile(baseDir,'neuronDataStruct.mat'),'neuronDataStruct');
neuronDataStruct = S.neuronDataStruct;
tsSec = cellfun(@(x) double(x)/fsNeur, {neuronDataStruct.timeStamps}, 'UniformOutput', false);

%% ---- load behavior labels + sync mapping ----
U = load(fullfile(baseDir,'UMAP.mat'));
V = load(fullfile(baseDir,'VideoSyncFrames.mat'));

% origDownsampEMGInd (required to map labels to emg indices)
if isfield(U,'origDownsampEMGInd')
    origDownsampEMGInd = U.origDownsampEMGInd(:);
else
    error('UMAP.mat is missing origDownsampEMGInd');
end

% choose label vector
switch lower(labelType)
    case "umap"
        if isfield(U,'regionAssignmentsFiltered')
            labels = U.regionAssignmentsFiltered;
        elseif isfield(U,'regionAssignmentFiltered')
            labels = U.regionAssignmentFiltered;
        else
            error('UMAP.mat is missing regionAssignment(s)Filtered');
        end
    otherwise
        error('labelType "%s" not supported yet. use "umap".', labelType);
end

% sync variables
if isfield(V,'frameNeuropixelSamples')
    frameNeuropixelSamples = V.frameNeuropixelSamples;
else
    error('VideoSyncFrames.mat is missing frameNeuropixelSamples');
end
if isfield(V,'frameEMGSamples')
    frameEMGSamples = V.frameEMGSamples;
else
    error('VideoSyncFrames.mat is missing frameEMGSamples');
end

% map emg indices -> neuropixels time base (ms units in your original code)
emgNeurSlope = (round(frameNeuropixelSamples{1}{end}(end)/30) - round(frameNeuropixelSamples{1}{1}(1)/30)) / ...
               (round(frameEMGSamples{1}{end}(end)/20)     - round(frameEMGSamples{1}{1}(1)/20));
emgNeurOffset = round(frameNeuropixelSamples{1}{1}(1)/30) - emgNeurSlope*round(frameEMGSamples{1}{1}(1)/20);

neurInds_ms = origDownsampEMGInd * emgNeurSlope + emgNeurOffset;  % neural ms indices
labelTimes_s = double(neurInds_ms) / 1000;                        % seconds

%% ---- make figure ----
figure('Color','w','Position',[100 100 1000 500]);
ax = axes; hold on;

% 1) behavior patches first (behind raster)
yl = [0, numel(tsSec)+1];
regs = unique(labels(~isnan(labels)));
patchColors = jet(numel(regs));

for r = 1:numel(regs)
    regVal = regs(r);
    idx = find(labels == regVal);

    if isempty(idx), continue; end

    starts = [1; find(diff(idx)~=1)+1];
    stops  = [find(diff(idx)~=1); numel(idx)];

    for b = 1:numel(starts)
        tStart = labelTimes_s(idx(starts(b)));
        tStop  = labelTimes_s(idx(stops(b)));

        % only draw if overlaps plotting window
        if tStop < t0 || tStart > t1
            continue;
        end

        % clip to window
        tStart = max(tStart, t0);
        tStop  = min(tStop,  t1);

        if tStop > tStart
            patch([tStart tStart tStop tStop], [yl(1) yl(2) yl(2) yl(1)], patchColors(r,:), ...
                'EdgeColor','none', 'FaceAlpha',0.15);
        end
    end
end

% 2) raster on top
for iNeuron = 1:numel(tsSec)
    t = tsSec{iNeuron};
    t = t(t>=t0 & t<=t1);
    if ~isempty(t)
        plot(t, iNeuron*ones(size(t)), 'k.', 'MarkerSize', 2);
    end
end

xlim([t0 t1]);
ylim(yl);
xlabel('Time (s)');
ylabel('Neuron index');
title('Raster snippet with behavior labels (UMAP patches)');
box off;

end
