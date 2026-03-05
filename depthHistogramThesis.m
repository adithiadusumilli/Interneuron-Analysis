function depthHistogramThesis(baseDir)
% depthhistogramthesis

% makes a cortex-only depth histogram split by pyramidal vs interneuron
% (mirrored bars: pyramidal right, interneuron left)

% input: basedir: session ProcessedDdata folder

% run j this: depthHistogramThesis("Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData")

arguments
    baseDir (1,1) string
end

%% ---- settings ----
conslidatedDataFoler = "X:\David\AnalysesData";

% define known base folders (match order used when AA_classifications.mat was created)
animalFolders = {
    'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData', ...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData', ...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData', ...
    'X:\David\ArenaRecordings\D043-020425-ArenaRecording\ProcessedData'
};

%% ---- load neuron depths ----
S = load(fullfile(baseDir,'neuronDataStruct.mat'),'neuronDataStruct');
neuronDataStruct = S.neuronDataStruct;

if ~isfield(neuronDataStruct,'depth')
    error('neuronDataStruct is missing the field "depth".');
end

depthAll = double([neuronDataStruct.depth]); % 1 x nUnits
nUnits = numel(depthAll);

%% ---- load cortex inds (you said cortexInds lives here) ----
F = load(fullfile(baseDir,'NeuralFiringRates1msBins10msGauss.mat'),'cortexInds');
if ~isfield(F,'cortexInds') || isempty(F.cortexInds)
    error('missing cortexInds in NeuralFiringRates1msBins10msGauss.mat');
end

cortexInds = double(F.cortexInds(:)');  % indices into neuronDataStruct
cortexInds = cortexInds(cortexInds >= 1 & cortexInds <= nUnits);

if isempty(cortexInds)
    error('cortexInds was empty after bounds-checking against neuronDataStruct.');
end

%% ---- load classifications (4x2 cell) ----
C = load(fullfile(conslidatedDataFoler,'AA_classifications.mat'),'classifications');
if ~isfield(C,'classifications') || isempty(C.classifications)
    error('AA_classifications.mat missing "classifications".');
end
classifications = C.classifications;

%% ---- match baseDir to the row used when the file was created ----
matchRow = find(contains(animalFolders, char(baseDir)), 1);

if isempty(matchRow)
    error('could not match baseDir to animalFolders. check the list in this function.');
end

% extract classification for this animal only
thisCortexLabels = classifications{matchRow, 1}; % 0=pyr, 1=int, nan=other

if isempty(thisCortexLabels) || numel(thisCortexLabels) ~= nUnits
    error('cortex classification vector size mismatch. expected 1x%d, got 1x%d.', ...
        nUnits, numel(thisCortexLabels));
end

%% ---- cortex-only pyramidal vs interneuron masks ----
isCortex = false(1, nUnits);
isCortex(cortexInds) = true;

pyrMask = isCortex & (thisCortexLabels == 0);
intMask = isCortex & (thisCortexLabels == 1);

depthPyr = depthAll(pyrMask);
depthInt = depthAll(intMask);

%% ---- plot mirrored depth histogram ----
binEdges = 0:100:4200; % 100 um bins; tweak if you want
binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;

nPyr = histcounts(depthPyr, binEdges);
nInt = histcounts(depthInt, binEdges);

figure('Color','w','Position',[100 100 520 520]); hold on;

% pyramidal to the right, interneurons to the left
barh(binCenters, nPyr, 1.0, 'EdgeColor','none');
barh(binCenters, -nInt, 1.0, 'EdgeColor','none');

set(gca,'YDir','reverse'); % depth increases downward like probe plots

xlabel('# of units');
ylabel('depth (\mum)');
title(sprintf('cortex depth distribution | pyr vs int | n=%d (pyr=%d, int=%d)', ...
    numel(cortexInds), sum(pyrMask), sum(intMask)));

% show positive tick labels on both sides
xmax = max([nPyr(:); nInt(:); 1]);
xlim([-xmax xmax]);
xt = unique(round(linspace(-xmax, xmax, 5)));
xticks(xt);
xticklabels(string(abs(xt)));

legend({'pyramidal','interneuron'}, 'Location','southeast');
box off;

end
