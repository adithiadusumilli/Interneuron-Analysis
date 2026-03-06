function depthHistogramThesis(baseDir)
% for fig 2!
% plotcortexdepthhistogrampyrint

% makes a cortex-only depth histogram using neuronDataStruct.depth and saved gmm spike-width classifications (0=pyramidal, 1=interneuron)


% j run depthHistogramThesis("Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData")

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


%% ---- load neuron depths ----
S = load(fullfile(baseDir,'neuronDataStruct.mat'),'neuronDataStruct');
neuronDataStruct = S.neuronDataStruct;

if ~isfield(neuronDataStruct,'depth')
    error('neuronDataStruct does not contain a depth field.');
end

allDepths_raw = double([neuronDataStruct.depth]); % 1 x nUnits (largest = brain surface)
nUnits = numel(allDepths_raw);

% convert so 0 µm corresponds to the surface and depth increases downward
allDepths = max(allDepths_raw) - allDepths_raw;

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

% depth binning (um)
binSize = 200;
maxDepth = ceil(max(cortexDepths)/binSize)*binSize;
binEdges = 0:binSize:maxDepth;

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
b = barh(ax, Y, X, 'stacked');

% set colors
b(1).FaceColor = [0.2 0.4 0.8]; % pyr
b(2).FaceColor = [0.9 0.4 0.2]; % int
b(1).EdgeColor = 'none';
b(2).EdgeColor = 'none';
b(1).FaceAlpha = 0.9;
b(2).FaceAlpha = 0.9;

% formatting to match example style
%set(ax,'YDir','reverse'); % 0 at top, increasing down
xlabel(ax,'# of units');
ylabel(ax,'Depth (\mum)');

title(ax, sprintf('Cortex Depth Distribution | Pyr vs Int | n=%d (pyr=%d, int=%d)', nTot, nPyr, nInt), 'FontSize',20);

legend(ax, {'Pyramidal','Interneuron'}, 'Location','southeast');
box(ax,'off');
ax.TickDir = 'out';
ax.LineWidth = 1;
ax.FontSize = 18;

end
