function paperDepthHistogram(baseDir)
% Cortex-only stacked depth histogram (0=pyramidal, 1=interneuron).
% Supported: D026, D020, D024, D050, D054. D043 is excluded.
% Uses getMouseDataNames when available to locate session data, including
% region-specific CFA firing-rate files for D050 and D054.
%
% Example:
% depthHistogramThesis("X:\David\ArenaRecordings\D050-121825-ArenaRecording\ProcessedData")

arguments
    baseDir (1,1) string
end

%% ---- settings ----
classificationFile = "X:\David\AnalysesData\AA_classifications.mat";
saveDir = "X:\David\AnalysesData\InterneuronAnalyses\AA Paper Plots\classifications";
fontSize = 10;
axesLineWidth = 0.5;
tickLength = [0.025 0.025];
numBins = 10;

% Classification rows in the saved six-session dataset, if the MAT file
% does not contain explicit mouseIDs/session metadata.
classificationOrder = ["D026","D020","D024","D043","D050","D054"];

%% ---- identify animal ----
match = regexp(char(baseDir), 'D\d{3}', 'match', 'once');
if isempty(match)
    error('Could not determine mouse ID from baseDir: %s', baseDir);
end
mouseID = string(match);
if mouseID == "D043"
    fprintf('D043 excluded: no histogram created.\n');
    return
end
if ~ismember(mouseID, classificationOrder)
    error('Unsupported mouse ID: %s', mouseID);
end

sessionName = regexp(char(baseDir), 'D\d{3}-\d{6}-ArenaRecording', 'match', 'once');
if isempty(sessionName)
    error('Could not determine session name from baseDir: %s', baseDir);
end
probeRegion = 'CFA';

%% ---- resolve files using David''s helper where available ----
% We do not assume a particular getMouseDataNames return-field spelling.
% Its returned paths are used when identifiable; otherwise the session
% folder and known file names are used.
neuronFile = fullfile(baseDir, 'neuronDataStruct.mat');
frFile = fullfile(baseDir, 'NeuralFiringRates1msBins10msGauss.mat');
if ismember(mouseID, ["D050","D054"])
    frFile = fullfile(baseDir, 'NeuralFiringRates1msBins10msGauss_CFA.mat');
end

if exist('getMouseDataNames', 'file') == 2
    try
        names = getMouseDataNames(char(mouseID), sessionName, probeRegion);
        if isstruct(names)
            fields = fieldnames(names);
            for k = 1:numel(fields)
                val = names.(fields{k});
                if ~(ischar(val) || (isstring(val) && isscalar(val)))
                    continue
                end
                val = string(val);
                if ~isfile(val)
                    continue
                end
                [~,nm,ext] = fileparts(val);
                if strcmpi(nm+ext, 'neuronDataStruct.mat')
                    neuronFile = val;
                elseif contains(lower(nm), 'neuralfiringrates1msbins10msgauss')
                    if ~ismember(mouseID, ["D050","D054"]) || contains(upper(nm),'CFA')
                        frFile = val;
                    end
                end
            end
        end
    catch ME
        warning('getMouseDataNames lookup failed (%s). Using session paths.', ME.message);
    end
end

if ~isfile(neuronFile)
    error('Cannot find neuronDataStruct file: %s', neuronFile);
end
if ~isfile(frFile)
    error('Cannot find firing-rate file: %s', frFile);
end

%% ---- load depths and cortex indices ----
S = load(neuronFile, 'neuronDataStruct');
if ~isfield(S, 'neuronDataStruct') || ~isfield(S.neuronDataStruct, 'depth')
    error('Missing neuronDataStruct.depth in %s', neuronFile);
end
rawDepths = double([S.neuronDataStruct.depth]);
nUnits = numel(rawDepths);
allDepths = max(rawDepths) - rawDepths;

F = load(frFile, 'cortexInds');
if ~isfield(F, 'cortexInds') || isempty(F.cortexInds)
    error('Missing cortexInds in %s', frFile);
end
cortexInds = double(F.cortexInds(:))';
if any(~isfinite(cortexInds) | cortexInds ~= round(cortexInds) | ...
       cortexInds < 1 | cortexInds > nUnits)
    error('Invalid cortexInds relative to neuronDataStruct (%d units).', nUnits);
end

%% ---- load GMM classifications and identify matching row ----
C = load(classificationFile);
if ~isfield(C,'classifications') || isempty(C.classifications)
    error('Missing classifications in %s', classificationFile);
end
classifications = C.classifications;

row = [];
for fieldName = ["mouseIDs","animalIDs"]
    if isfield(C, char(fieldName))
        ids = string(C.(char(fieldName)));
        row = find(strcmpi(ids(:), mouseID), 1);
        if ~isempty(row), break; end
    end
end
if isempty(row) && isfield(C,'baseSessionNames')
    ids = string(C.baseSessionNames);
    row = find(contains(ids(:), mouseID), 1);
end
if isempty(row)
    if size(classifications,1) ~= numel(classificationOrder)
        error(['Cannot safely match classification row: no mouseID metadata ' ...
            'and classifications has %d rows (expected %d).'], ...
            size(classifications,1), numel(classificationOrder));
    end
    row = find(classificationOrder == mouseID, 1);
    fprintf('Using six-session classification row order: %s -> row %d\n',mouseID,row);
end
if row > size(classifications,1)
    error('Classification row %d does not exist.',row);
end
labels = double(classifications{row,1});
labels = labels(:)';
if max(cortexInds) > numel(labels)
    error('cortexInds exceed length of cortex classifications (%d).',numel(labels));
end
cortexLabels = labels(cortexInds);
cortexDepths = allDepths(cortexInds);

%% ---- filter and bin ----
keep = isfinite(cortexDepths) & ismember(cortexLabels,[0 1]);
cortexDepths = cortexDepths(keep);
cortexLabels = cortexLabels(keep);
pyrDepths = cortexDepths(cortexLabels == 0);
intDepths = cortexDepths(cortexLabels == 1);
nPyr = numel(pyrDepths);
nInt = numel(intDepths);
nTot = nPyr + nInt;
if nTot == 0
    error('No cortex units with valid depth and pyramidal/interneuron labels.');
end
maxDepth = max(cortexDepths);
if maxDepth <= 0
    error('Cortex depth range is zero; cannot construct depth bins.');
end
binEdges = linspace(0,maxDepth,numBins+1);
binCenters = binEdges(1:end-1) + diff(binEdges)/2;
pyrCounts = histcounts(pyrDepths,binEdges);
intCounts = histcounts(intDepths,binEdges);

%% ---- plot ----
fig = figure('Name',sprintf('%s cortex depth histogram',mouseID), ...
    'Color','w','Position',[100 100 650 450]);
ax = axes(fig); hold(ax,'on');
b = barh(ax,binCenters(:),[pyrCounts(:),intCounts(:)],'stacked');
b(1).FaceColor = [0 0 1];
b(2).FaceColor = [1 0 0];
for k = 1:2
    b(k).EdgeColor = 'none';
    b(k).FaceAlpha = 0.75;
end
xlabel(ax,'Number of units','FontSize',fontSize);
ylabel(ax,'Depth (\mum)','FontSize',fontSize);
title(ax,sprintf('%s cortex depth distribution (n=%d; pyr=%d, int=%d)', ...
    mouseID,nTot,nPyr,nInt),'FontSize',fontSize);
legend(ax,{'Pyramidal neuron','Interneuron'},'Location','northeast', ...
    'FontSize',fontSize,'Box','off');
set(ax,'FontSize',fontSize,'LineWidth',axesLineWidth, ...
    'TickDir','out','TickLength',tickLength,'XColor','k','YColor','k', ...
    'Color','w');
grid(ax,'off'); box(ax,'off');

%% ---- save ----
if ~isfolder(saveDir), mkdir(saveDir); end
pdfFile = fullfile(saveDir,sprintf('%s_cortex_depth_histogram.pdf',mouseID));
figFile = fullfile(saveDir,sprintf('%s_cortex_depth_histogram.fig',mouseID));
exportgraphics(fig,pdfFile,'ContentType','vector','BackgroundColor','white');
savefig(fig,figFile);
fprintf('%s: %d pyramidal, %d interneurons.\nSaved: %s\n', ...
    mouseID,nPyr,nInt,pdfFile);
end
