function plotCortexDepthHistogramPyrInt(baseDir)
% plotcortexdepthhistogrampyrint

% makes a cortex-only depth histogram using neuronDataStruct.depth and saved gmm spike-width classifications (0=pyramidal, 1=interneuron)


% j run plotCortexDepthHistogramPyrInt("Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData")

arguments
    baseDir (1,1) string
end

%% ---- settings ----
conslidatedDataFoler = 'X:\David\AnalysesData';

% this order must match how AA_classifications.mat was created
animalFolders = {
    'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData', ...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData', ...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData', ...
    'X:\David\ArenaRecordings\D043-020425-ArenaRecording\ProcessedData'
};

% depth binning (um)
binSize = 200;
maxDepth = 4200;
binEdges = 0:binSize:maxDepth;

%% ---- load neuron depths ----
S = load(fullfile(baseDir,'neuronDataStruct.mat'),'neuronDataStruct');
neuronDataStruct = S.neuronDataStruct;

if ~isfield(neuronDataStruct,'depth')
    error('neuronDataStruct does not contain a depth field.');
end

allDepths = double([neuronDataStruct.depth]); % 1 x nUnits
nUnits = numel(allDepths);

%% ---- load cortex indices (these are indices into neuronDataStruct) ----
F = load(fullfile(baseDir,'NeuralFiringRates1msBins10msGauss.mat'),'cortexInds');
if ~isfield(F,'cortexInds') || isempty(F.cortexInds)
    error('could not find cortexInds in NeuralFiringRates1msBins10msGauss.mat');
end
cortexInds = double(F.cortexInds(:))';

% bounds check
cortexInds = cortexInds(cortexInds >= 1 & cortexInds <= nUnits);
if isempty(cortexInds)
    error('cortexInds was empty after bounds checking against neuronDataStruct (%d units).', nUnits);
end

%% ---- load classifications + match this animal ----
C = load(fullfile(conslidatedDataFoler,'AA_classifications.mat'),'classifications');
if ~isfield(C,'classifications') || isempty(C.classifications)
    error('AA_classifications.mat does not contain classifications.');
end
classifications = C.classifications;

matchRow = find(contains(string(animalFolders), baseDir), 1);
if isempty(matchRow)
    error('could not match baseDir to an entry in animalFolders. update animalFolders to include this session.');
end

% cortex labels are stored in column 1
cortexLabelsAll = classifications{matchRow,1};
if isempty(cortexLabelsAll)
    error('classifications{%d,1} was empty (cortex labels).', matchRow);
end
cortexLabelsAll = double(cortexLabelsAll);

% restrict to cortex units only
cortexLabels = cortexLabelsAll(cortexInds); % 0=pyr, 1=int, nan=other
cortexDepths = allDepths(cortexInds);

% keep only labeled pyr/int and non-nan depths
keep = ~isnan(cortexDepths) & ~isnan(cortexLabels) & (cortexLabels == 0 | cortexLabels == 1);
cortexDepths = cortexDepths(keep);
cortexLabels = cortexLabels(keep);

pyrDepths = cortexDepths(cortexLabels == 0);
intDepths = cortexDepths(cortexLabels == 1);

nPyr = numel(pyrDepths);
nInt = numel(intDepths);
nTot = nPyr + nInt;

if nTot == 0
    error('no cortex units with valid depth + pyr/int labels were found.');
end

%% ---- compute histogram counts per depth bin ----
pyrCounts = histcounts(pyrDepths, binEdges);
intCounts = histcounts(intDepths, binEdges);

binCenters = binEdges(1:end-1) + diff(binEdges)/2;

%% ---- plot (single-sided, stacked) ----
figure('Color','w','Position',[100 100 520 520]);
ax = axes; hold(ax,'on');

Y = binCenters(:);
X = [pyrCounts(:), intCounts(:)];

% stacked horizontal bars so everything stays on the right
barh(ax, Y, X, 'stacked');

% formatting to match example style
set(ax,'YDir','reverse'); % 0 at top, increasing down
xlabel(ax,'# of units');
ylabel(ax,'depth (\mum)');

title(ax, sprintf('cortex depth distribution | pyr vs int | n=%d (pyr=%d, int=%d)', nTot, nPyr, nInt));

legend(ax, {'pyramidal','interneuron'}, 'Location','southeast');
box(ax,'off');
ax.TickDir = 'out';
ax.LineWidth = 1;
ax.FontSize = 12;

end
