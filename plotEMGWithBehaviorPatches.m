function plotEMGWithBehaviorPatches(baseDir, t0, t1, emgChannel, labelType)
% plotEMGWithBehaviorPatches

% plots an emg snippet with behavior labels as background patches

% inputs:
%   baseDir     : session ProcessedData folder
%   t0, t1      : window in seconds
%   emgChannel  : which channel to plot (if downsampEMG is multichannel). use 1 if unsure.
%   labelType   : 'umap' (default)

% example:
%   plotEMGWithBehaviorPatches('Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData',600,675,1,'umap')

arguments
    baseDir (1,1) string
    t0 (1,1) double
    t1 (1,1) double
    emgChannel (1,1) double = 1
    labelType (1,1) string = "umap"
end

fsEmg = 1000;  % bc downsamp from 20khz

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
    emg = emgAll(:, emgChannel);
end

emgTime = (0:numel(emg)-1) / fsEmg;

%% ---- load behavior + mapping (same as raster function) ----
U = load(fullfile(baseDir,'UMAP.mat'));
V = load(fullfile(baseDir,'VideoSyncFrames.mat'));

if isfield(U,'origDownsampEMGInd')
    origDownsampEMGInd = U.origDownsampEMGInd(:);
else
    error('UMAP.mat is missing origDownsampEMGInd');
end

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

if ~isfield(V,'frameNeuropixelSamples') || ~isfield(V,'frameEMGSamples')
    error('VideoSyncFrames.mat missing frameNeuropixelSamples and/or frameEMGSamples');
end

frameNeuropixelSamples = V.frameNeuropixelSamples;
frameEMGSamples       = V.frameEMGSamples;

% map emg indices -> neuropixels ms indices, then to seconds (same mapping)
emgNeurSlope = (round(frameNeuropixelSamples{1}{end}(end)/30) - round(frameNeuropixelSamples{1}{1}(1)/30)) / ...
               (round(frameEMGSamples{1}{end}(end)/20)     - round(frameEMGSamples{1}{1}(1)/20));
emgNeurOffset = round(frameNeuropixelSamples{1}{1}(1)/30) - emgNeurSlope*round(frameEMGSamples{1}{1}(1)/20);

neurInds_ms = origDownsampEMGInd * emgNeurSlope + emgNeurOffset;
labelTimes_s = double(neurInds_ms) / 1000;

%% ---- plot ----
figure('Color','w','Position',[100 100 1000 350]);
ax = axes; hold on;

% behavior patches
yl = [-inf inf]; % will reset after plotting emg
regs = unique(labels(~isnan(labels)));
patchColors = jet(numel(regs));

% plot emg first to get y-limits
keep = emgTime>=t0 & emgTime<=t1;
plot(emgTime(keep), emg(keep), 'k');
xlim([t0 t1]);
yl = get(gca,'YLim');

for r = 1:numel(regs)
    regVal = regs(r);
    idx = find(labels == regVal);

    if isempty(idx), continue; end

    starts = [1; find(diff(idx)~=1)+1];
    stops  = [find(diff(idx)~=1); numel(idx)];

    for b = 1:numel(starts)
        tStart = labelTimes_s(idx(starts(b)));
        tStop  = labelTimes_s(idx(stops(b)));

        if tStop < t0 || tStart > t1
            continue;
        end

        tStart = max(tStart, t0);
        tStop  = min(tStop,  t1);

        if tStop > tStart
            patch([tStart tStart tStop tStop], [yl(1) yl(2) yl(2) yl(1)], patchColors(r,:), ...
                'EdgeColor','none','FaceAlpha',0.15);
        end
    end
end

% bring emg to front
uistack(findobj(gca,'Type','line'),'top');

xlabel('Time (s)');
ylabel('EMG (a.u.)');
title('EMG snippet with behavior labels (UMAP patches)');
box off;

end
