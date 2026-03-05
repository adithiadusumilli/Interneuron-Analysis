% example: depth histogram for one animal (cortex only), split by pyr vs int
% - pulls depth from neuronDataStruct.depth (um)
% - uses cortexInds from NeuralFiringRates1msBins10msGauss.mat
% - uses pyr/int labels from AA_classifications.mat (in consolidated folder)

% D024 path
  
baseDir = "Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData";
conslidatedDataFoler = "X:\David\AnalysesData"; % (keeping your variable name spelling)

%% ---- load depths ----
S = load(fullfile(baseDir,"neuronDataStruct.mat"), "neuronDataStruct");
neuronDataStruct = S.neuronDataStruct;

if ~isfield(neuronDataStruct, "depth")
    error('neuronDataStruct is missing the field "depth".');
end

depthAll = double([neuronDataStruct.depth]); % 1 x nUnits

%% ---- load cortex indices (this is where you said cortexInds lives) ----
F = load(fullfile(baseDir,"NeuralFiringRates1msBins10msGauss.mat"), "cortexInds");
if ~isfield(F,"cortexInds") || isempty(F.cortexInds)
    error('missing cortexInds in NeuralFiringRates1msBins10msGauss.mat');
end
cortexInds = double(F.cortexInds(:)');  % 1 x nCortex

% bounds check
nUnits = numel(depthAll);
cortexInds = cortexInds(cortexInds>=1 & cortexInds<=nUnits);

%% ---- load classifications (pyr vs int) ----
C = load(fullfile(conslidatedDataFoler,"AA_classifications.mat"));
% try common fields; edit here if your file uses different names
[pyrMaskAll, intMaskAll] = local_get_pyr_int_masks(C, nUnits);

% cortex-only masks
isCortex = false(1,nUnits);
isCortex(cortexInds) = true;

pyrCortex = isCortex & pyrMaskAll;
intCortex = isCortex & intMaskAll;

depthPyr = depthAll(pyrCortex);
depthInt = depthAll(intCortex);

%% ---- plot histogram like your figure (horizontal bars) ----
binEdges = 0:100:4200; % 100 um bins; tweak if you want
binCenters = (binEdges(1:end-1)+binEdges(2:end))/2;

nPyr = histcounts(depthPyr, binEdges);
nInt = histcounts(depthInt, binEdges);

figure('Color','w','Position',[100 100 520 520]); hold on;

% plot pyramidal to the right, interneurons to the left (mirrored)
barh(binCenters,  nPyr, 1.0, 'FaceAlpha',0.85, 'EdgeColor','none');     % pyr
barh(binCenters, -nInt, 1.0, 'FaceAlpha',0.85, 'EdgeColor','none');     % int

set(gca,'YDir','reverse'); % depth increases downward like probe plots

xlabel('# of units');
ylabel('depth (\mum)');
title('cortex unit depth distribution | pyramidal vs interneuron');

% make x-axis show positive counts on both sides
xmax = max([nPyr(:); nInt(:); 1]);
xlim([-xmax xmax]);
xticks(-xmax:round(xmax/4):xmax);
xticklabels(string(abs(xticks)));

legend({'pyramidal','interneuron'}, 'Location','southeast');
box off;

%% ---------------- local helper ----------------
function [pyrMask, intMask] = local_get_pyr_int_masks(C, nUnits)
% returns logical masks over ALL units (length nUnits)

% initialize
pyrMask = false(1,nUnits);
intMask = false(1,nUnits);

% candidate field names people commonly use
pyrFields = {'pyrInds','pyr_inds','pyrIdx','pyrUnits','pyramidalInds','wideInds','wideSpikingInds'};
intFields = {'intInds','int_inds','intIdx','intUnits','interneuronInds','narrowInds','narrowSpikingInds'};

% helper to grab a field if it exists (also if nested under "classifications")
valPyr = local_find_field(C, pyrFields);
valInt = local_find_field(C, intFields);

if isempty(valPyr) || isempty(valInt)
    error(['could not find pyramidal/interneuron fields in AA_classifications.mat. ' ...
           'open it with: whos("-file",fullfile(conslidatedDataFoler,"AA_classifications.mat")) ' ...
           'and tell me the field names you see for pyr/int.']);
end

pyrMask = local_as_mask(valPyr, nUnits);
intMask = local_as_mask(valInt, nUnits);

% if they overlap, prefer explicit int/pyr but warn
if any(pyrMask & intMask)
    warning('pyr and int masks overlap. check classification fields.');
end
end

function val = local_find_field(S, candidates)
val = [];
for i = 1:numel(candidates)
    f = candidates{i};
    if isfield(S,f) && ~isempty(S.(f))
        val = S.(f);
        return;
    end
end
% nested
if isfield(S,'classifications') && isstruct(S.classifications)
    for i = 1:numel(candidates)
        f = candidates{i};
        if isfield(S.classifications,f) && ~isempty(S.classifications.(f))
            val = S.classifications.(f);
            return;
        end
    end
end
end

function mask = local_as_mask(v, nUnits)
if islogical(v)
    mask = v(:)';
    if numel(mask) ~= nUnits
        error('logical mask length (%d) does not match nUnits (%d).', numel(mask), nUnits);
    end
else
    inds = double(v(:)');
    inds = inds(inds>=1 & inds<=nUnits);
    mask = false(1,nUnits);
    mask(inds) = true;
end
end
