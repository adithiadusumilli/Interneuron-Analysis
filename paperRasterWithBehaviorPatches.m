function paperRasterWithBehaviorPatches(baseDir, t0, t1, labelType)
% Plots a cortex-only raster snippet with behavior labels as background patches.
% Interneuron spikes = red; pyramidal neuron spikes = black.
% Behavior patches retain their original 30%-opacity appearance.
% Automatically saves a vector PDF.

% RUN:
% paperRasterWithBehaviorPatches("Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData", 600, 675, "classifier")

arguments
    baseDir (1,1) string
    t0 (1,1) double
    t1 (1,1) double
    labelType (1,1) string = "classifier"
end

fsNeur = 30000;

if t1 <= t0
    error('t1 must be > t0');
end

%% ---- load cortex indices from firing rate file ----

F = load(fullfile(baseDir, ...
    'NeuralFiringRates1msBins10msGauss.mat'), 'cortexInds');

if ~isfield(F,'cortexInds') || isempty(F.cortexInds)
    error('Could not find cortexInds in firing rate file.');
end

cortexInds = F.cortexInds;

%% ---- load spikes ----

S = load(fullfile(baseDir,'neuronDataStruct.mat'), ...
    'neuronDataStruct');

neuronDataStruct = S.neuronDataStruct;

tsSec_all = cellfun( ...
    @(x) double(x)/fsNeur, ...
    {neuronDataStruct.timeStamps}, ...
    'UniformOutput',false);

nAll = numel(tsSec_all);

if islogical(cortexInds)
    cortexInds = find(cortexInds);
end

cortexInds = cortexInds(:)';
cortexInds = cortexInds( ...
    cortexInds >= 1 & cortexInds <= nAll);

if isempty(cortexInds)
    error('cortexInds is empty after bounds checking.');
end

tsSec = tsSec_all(cortexInds);
nNeur = numel(tsSec);

%% ---- load saved neuron classifications ----

classificationFile = ...
    "X:\David\AnalysesData\AA_classifications.mat";

if ~isfile(classificationFile)
    error('Classification file not found: %s',classificationFile);
end

C = load(classificationFile,'mouseIDs','classifications');

% Extract mouse ID from session path.
mouseMatch = regexp(char(baseDir),'D\d+','match','once');

if isempty(mouseMatch)
    error('Could not determine mouse ID from baseDir.');
end

mouseID = string(mouseMatch);

mouseRow = find(string(C.mouseIDs) == mouseID,1,'first');

if isempty(mouseRow)
    error('Mouse %s not found in AA_classifications.mat.',mouseID);
end

% Column 1 contains cortex classifications.
% 0 = pyramidal neuron
% 1 = interneuron
% NaN = unclassified

% Saved classification labels may be indexed by all recorded neurons.
% Match the cortex raster ordering using cortexInds when necessary.
allNeuronTypes = double(C.classifications{mouseRow,1}(:));

if numel(allNeuronTypes) == nAll
    neuronTypes = allNeuronTypes(cortexInds);
elseif numel(allNeuronTypes) == nNeur
    neuronTypes = allNeuronTypes;
else
    error(['Classification count mismatch for %s: ' ...
        '%d labels, %d total neurons, %d cortex neurons.'], ...
        mouseID,numel(allNeuronTypes),nAll,nNeur);
end
neuronTypes = neuronTypes(:);

% Reassign raster y-axis positions by neuron type without changing spike data.
% Interneurons (red) at bottom, pyramidal neurons (black) above.
% Unclassified units (if any) are placed at the top.
interneuronIdx = find(neuronTypes == 1);
pyramidalIdx = find(neuronTypes == 0);
unclassifiedIdx = find(~ismember(neuronTypes,[0 1]));
newOrder = [interneuronIdx; pyramidalIdx; unclassifiedIdx];
rasterY = zeros(nNeur,1);
rasterY(newOrder) = 1:nNeur;

fprintf('\n%s raster classifications:\n',mouseID);
fprintf('Pyramidal neurons: %d\n',sum(neuronTypes == 0));
fprintf('Interneurons: %d\n',sum(neuronTypes == 1));
fprintf('Unclassified: %d\n',sum(isnan(neuronTypes)));

%% ---- load behavior labels EXACTLY as in the non-chunked analysis ----
% runNonChunkedXcorrPerBehavior_saveOnly_fullperm uses the original
% origDownsampEMGInd directly as 1-ms firing-rate bin indices.
% No VideoSyncFrames transformation is applied.

U = load(fullfile(baseDir,'UMAP.mat'), ...
    'origDownsampEMGInd','classifierLabels','classifierBehvs');

required = {'origDownsampEMGInd','classifierLabels','classifierBehvs'};
for k = 1:numel(required)
    if ~isfield(U,required{k}) || isempty(U.(required{k}))
        error('UMAP.mat missing/empty field: %s',required{k});
    end
end

if lower(labelType) ~= "classifier"
    error(['This version matches the saved per-behavior analysis, ' ...
        'which only supports labelType="classifier".']);
end

manBehvNames = {'climbdown','climbup','eating','grooming', ...
    'jumpdown','jumping','rearing','still','walkflat','walkgrid'};
nColors = 11;

canonicalLookup = containers.Map;
canonicalLookup('climbdown') = 1;
canonicalLookup('climbup') = 2;
canonicalLookup('eating') = 3;
canonicalLookup('eat') = 3;
canonicalLookup('grooming') = 4;
canonicalLookup('groom') = 4;
canonicalLookup('jumpdown') = 5;
canonicalLookup('jumping') = 6;
canonicalLookup('jumpacross') = 6;
canonicalLookup('rearing') = 7;
canonicalLookup('rear') = 7;
canonicalLookup('still') = 8;
canonicalLookup('walkflat') = 9;
canonicalLookup('walkgrid') = 10;

classifierBehvs = cellstr(string(U.classifierBehvs(:)));
classBehvNumbers = zeros(1,numel(classifierBehvs));
for k = 1:numel(classifierBehvs)
    cleanName = lower(strrep(strrep(classifierBehvs{k},' ',''),'_',''));
    if isKey(canonicalLookup,cleanName)
        classBehvNumbers(k) = canonicalLookup(cleanName);
    else
        classBehvNumbers(k) = 0;
    end
end

origInd = double(U.origDownsampEMGInd(:));
rawLabels = double(U.classifierLabels(:));
n = min(numel(origInd),numel(rawLabels));
origInd = origInd(1:n);
rawLabels = rawLabels(1:n);
labelsReduced = nan(n,1);
for k = 1:n
    v = rawLabels(k);
    if isnan(v)
        continue;
    elseif v == 0
        labelsReduced(k) = 0;
    elseif v >= 1 && v <= numel(classBehvNumbers) && v == floor(v)
        labelsReduced(k) = classBehvNumbers(v);
    end
end

% Each index corresponds to a 1-ms neural bin. Use bin START time for
% raster patches, and the following millisecond as the bin END time.
valid = isfinite(origInd) & origInd == floor(origInd) & origInd >= 1 & ...
    isfinite(labelsReduced);
origInd = origInd(valid);
labelsReduced = labelsReduced(valid);

% Sort and remove repeated bin indices, as the analysis assigns one label
% per bin (last assignment wins).
[origInd,ord] = sort(origInd);
labelsReduced = labelsReduced(ord);
[origInd,keep] = unique(origInd,'last');
labelsReduced = labelsReduced(keep);

startTimes = (origInd-1)/1000;
endTimes = origInd/1000;
inWin = endTimes > t0 & startTimes < t1;
patchBinInd = origInd(inWin);
patchLabels = labelsReduced(inWin);

fprintf('Behavior bins overlapping [%.3f, %.3f] s: %d\n', ...
    t0,t1,numel(patchBinInd));
if ~isempty(patchBinInd)
    fprintf('Behavior codes in window: ');
    fprintf('%g ',unique(patchLabels));
    fprintf('\n');
else
    warning('No classifier behavior bins overlap this raster time window.');
end

%% ---- create figure ----

fig = figure( ...
    'Color','w', ...
    'Position',[100 100 1100 520]);

ax = axes(fig);
hold(ax,'on');

yl = [0,nNeur+1];

cmap = local_behavior_cmap(nColors);

%% ---- behavior background patches ----
% Draw original palette with original FaceAlpha = 0.30. Group only
% CONSECUTIVE 1-ms bins with the same label, so gaps are not filled.
if ~isempty(patchBinInd)
    blockStarts = find([true; diff(patchBinInd) ~= 1 | diff(patchLabels) ~= 0]);
    blockEnds = [blockStarts(2:end)-1; numel(patchBinInd)];

    for b = 1:numel(blockStarts)
        i0 = blockStarts(b);
        i1 = blockEnds(b);
        beh = patchLabels(i0);
        cInd = local_label_to_color_index(beh,nColors);
        if isnan(cInd)
            continue;
        end
        tStart = max((patchBinInd(i0)-1)/1000,t0);
        tStop = min(patchBinInd(i1)/1000,t1);
        if tStop > tStart
            patch(ax,[tStart tStart tStop tStop], ...
                [yl(1) yl(2) yl(2) yl(1)],cmap(cInd,:), ...
                'EdgeColor','none','FaceAlpha',0.30);
        end
    end
    fprintf('Plotted %d behavior blocks.\n',numel(blockStarts));
end

%% ---- raster spikes colored by neuron type ----

tickHalfHeight = 0.35;

interneuronColor = [1 0 0];  % red
pyramidalColor = [0 0 0];   % black

for iNeuron = 1:nNeur

    t = tsSec{iNeuron};
    t = t(t >= t0 & t <= t1);

    if isempty(t)
        continue;
    end

    % 1 = interneuron; 0 = pyramidal.
    % Unclassified units remain black.
    if neuronTypes(iNeuron) == 1
        spikeColor = interneuronColor;
    else
        spikeColor = pyramidalColor;
    end

    x = [t(:)';t(:)'];

    y = [ ...
        (rasterY(iNeuron)-tickHalfHeight)*ones(1,numel(t)); ...
        (rasterY(iNeuron)+tickHalfHeight)*ones(1,numel(t))];

    line(ax,x,y, ...
        'Color',spikeColor, ...
        'LineWidth',0.5);
end

%% ---- axes formatting (original style retained) ----

xlim(ax,[t0 t1]);
ylim(ax,yl);

ax.FontSize = 14;
ax.TickDir = 'out';
ax.LineWidth = 1;

xlabel(ax,'Time (s)','FontSize',16);
ylabel(ax,'M1 Neuron Index','FontSize',16);

title(ax, ...
    sprintf('M1 Raster with %s Behavior Patches | %.1f–%.1f s | n=%d units', ...
    labelType,t0,t1,nNeur), ...
    'FontSize',18);

box(ax,'off');

%% ---- behavior legend ----

behNames = local_behavior_names(U,manBehvNames,nColors);

h = gobjects(1,numel(behNames));

for k = 1:numel(behNames)

    h(k) = patch(ax, ...
        nan,nan,cmap(k,:), ...
        'EdgeColor','none', ...
        'FaceAlpha',0.30);
end

legend(ax,h,cellstr(behNames), ...
    'Location','eastoutside');

%% ---- save vector PDF ----

saveDir = ...
    "X:\David\AnalysesData\InterneuronAnalyses\AA Paper Plots\rasters";

if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

safeLabel = char(regexprep(labelType,'[^a-zA-Z0-9]','_'));

fileBase = sprintf( ...
    '%s_Raster_%s_%.1f_to_%.1f_s', ...
    char(mouseID),safeLabel,t0,t1);

pdfFile = fullfile(saveDir,[fileBase '.pdf']);

% Confirm that the original figure is still valid before exporting.
if ~isgraphics(fig,'figure')
    error('Raster figure was deleted before PDF export.');
end

drawnow;

% Export the axes, avoiding figure-level export issues with some MATLAB
% versions while preserving the behavior legend outside the axes.
try
    exportgraphics(fig,pdfFile, ...
        'ContentType','vector', ...
        'BackgroundColor','white');
catch ME
    if ~isgraphics(fig,'figure')
        error('Raster figure was deleted before PDF export: %s',ME.message);
    end
    warning('paperRaster:FigureExportFailed', ...
        'Figure export failed (%s). Retrying with axes.',ME.message);
    exportgraphics(ax,pdfFile, ...
        'ContentType','vector', ...
        'BackgroundColor','white');
end

fprintf('\nSaved vector PDF:\n%s\n',pdfFile);

end

%% ================================================================
%% LOCAL HELPERS
%% ================================================================

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
        0.00 0.45 0.74;
        0.85 0.33 0.10;
        0.00 0.62 0.45;
        0.93 0.69 0.13;
        0.49 0.18 0.56;
        0.30 0.75 0.93;
        1.00 0.50 0.00
    ];
end

end

function cInd = local_label_to_color_index(beh,nColors)

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

    cInd = beh+1;
end

end

function behNames = local_behavior_names(U,manBehvNames,nColors)

if nColors == 7

    if isfield(U,'regionBehvAssignments') && ...
            ~isempty(U.regionBehvAssignments)

        tmp = string(U.regionBehvAssignments(:)');
        tmp = tmp(1:min(7,numel(tmp)));

        if numel(tmp) < 7
            tmp(end+1:7) = ...
                "umap region " + (numel(tmp)+1:7);
        end

        behNames = tmp;

    else

        behNames = "umap region " + (1:7);
    end

else

    behNames = ["unlabeled",string(manBehvNames)];
end

end

function out = iff(cond,a,b)

if cond
    out = a;
else
    out = b;
end

end
