function mapTransitionsToCanonicalBehaviors(mouseID, baseSessionName, probeRegion)
% mapTransitionsToCanonicalBehaviors

% Uses David's getMouseDataNames() to locate the processed-data folder, then maps EMG transition windows to canonical behavior labels.

% For the newer animals, there are no separate manual/UMAP behavior labels in the old format, so the manual and region-label outputs are populated from the classifier labels to preserve downstream pipeline compatibility.

% INPUTS
%   mouseID 
%   baseSessionName session name used by getMouseDataNames
%   probeRegion aka probe-region argument used by getMouseDataNames

% EXAMPLE
%   mapTransitionsToCanonicalBehaviors('D050', baseSessionName, probeRegion)

% This function:
%   - uses getMouseDataNames to locate the animal/session data
%   - loads EMG transition indices for each muscle channel
%   - maps each transition index into the reduced behavior-label time axis
%   - assigns behavior labels using a consistent canonical scheme

% Canonical output:
%   regionLabelsPerTransition      0-10
%   manualLabelsPerTransition      0-10
%   classifierLabelsPerTransition  0-10

% Canonical behavior ordering:
%   1  climbdown
%   2  climbup
%   3  eating
%   4  grooming
%   5  jumpdown
%   6  jumping
%   7  rearing
%   8  still
%   9  walkflat
%   10 walkgrid

% Behaviors not represented in the canonical list are mapped to 0.

%% ============================================================
%% Locate files using David's getMouseDataNames
%% ============================================================

dataNames = getMouseDataNames(mouseID, baseSessionName, probeRegion);

if ~isfield(dataNames, 'processedDataFolder') || isempty(dataNames.processedDataFolder)

    error(['getMouseDataNames did not return a valid ' 'processedDataFolder for %s.'], mouseID);
end

folderPath = dataNames.processedDataFolder;

fprintf('\n========================================\n');
fprintf('mapping transition behaviors for %s\n', mouseID);
fprintf('========================================\n');
fprintf('processed data folder:\n%s\n\n', folderPath);

%% ============================================================
%% Files
%% ============================================================

emgFile = fullfile(folderPath, 'EMG_Neural_AllChannels.mat');
umapFile = fullfile(folderPath, 'UMAP.mat');

if ~isfile(emgFile)
    error('Missing EMG file:\n%s', emgFile);
end

if ~isfile(umapFile)
    error('Missing UMAP file:\n%s', umapFile);
end

%% ============================================================
%% Load transition indices
%% ============================================================

E = load(emgFile, 'validTransitionsCell');

if ~isfield(E, 'validTransitionsCell')
    error('%s does not contain validTransitionsCell.', emgFile);
end

validTransitionsCell = E.validTransitionsCell;

%% ============================================================
%% Load behavior variables
%% ============================================================

U = load(umapFile, ...
    'origDownsampEMGInd', ...
    'classifierLabels', ...
    'classifierBehvs', ...
    'regionBehvAssignments');

requiredFields = { ...
    'origDownsampEMGInd', ...
    'classifierLabels', ...
    'classifierBehvs'};

for iField = 1:numel(requiredFields)

    thisField = requiredFields{iField};

    if ~isfield(U, thisField)
        error('%s does not contain required variable %s.', ...
            umapFile, thisField);
    end
end

origDownsampEMGInd = U.origDownsampEMGInd;
classifierLabels = U.classifierLabels;
classifierBehvs = U.classifierBehvs;

% Preserve this variable for downstream compatibility if it exists.
if isfield(U, 'regionBehvAssignments')
    regionBehvAssignments = U.regionBehvAssignments;
else
    regionBehvAssignments = [];
end

%% ============================================================
%% Canonical behavior definitions
%% ============================================================

manBehvNames = { ...
    'climbdown', ...
    'climbup', ...
    'eating', ...
    'grooming', ...
    'jumpdown', ...
    'jumping', ...
    'rearing', ...
    'still', ...
    'walkflat', ...
    'walkgrid'};

%% ============================================================
%% Lookup table for alternate names
%% ============================================================

canonicalLookup = containers.Map;

canonicalLookup('climbdown')  = 1;
canonicalLookup('climbup')    = 2;

canonicalLookup('eating')     = 3;
canonicalLookup('eat')        = 3;

canonicalLookup('grooming')   = 4;
canonicalLookup('groom')      = 4;

canonicalLookup('jumpdown')   = 5;

canonicalLookup('jumping')    = 6;
canonicalLookup('jumpacross') = 6;

canonicalLookup('rearing')    = 7;
canonicalLookup('rear')       = 7;

canonicalLookup('still')      = 8;

canonicalLookup('walkflat')   = 9;
canonicalLookup('walkgrid')   = 10;

%% ============================================================
%% Build manual-label remap
%% ============================================================
% For the new animals, the "manual" labels are deliberately copied
% from the classifier labeling to preserve compatibility with older
% downstream analysis code.
%
% Per-animal classifier index -> canonical index 0:10.

nManualBehv = numel(classifierBehvs);

manBehvNumbers = zeros(1, nManualBehv);

for iBehv = 1:nManualBehv

    thisName = classifierBehvs{iBehv};

    cleanName = lower( ...
        strrep( ...
        strrep(thisName, ' ', ''), ...
        '_', ''));

    if isKey(canonicalLookup, cleanName)

        manBehvNumbers(iBehv) = ...
            canonicalLookup(cleanName);

    else

        % Behavior not represented in canonical list.
        manBehvNumbers(iBehv) = 0;

    end
end

%% ============================================================
%% Build classifier-label remap
%% ============================================================

nClassBehv = numel(classifierBehvs);

classBehvNumbers = zeros(1, nClassBehv);

for iBehv = 1:nClassBehv

    thisName = classifierBehvs{iBehv};

    cleanName = lower( ...
        strrep( ...
        strrep(thisName, ' ', ''), ...
        '_', ''));

    if isKey(canonicalLookup, cleanName)

        classBehvNumbers(iBehv) = ...
            canonicalLookup(cleanName);

    else

        classBehvNumbers(iBehv) = 0;

    end
end

%% ============================================================
%% Print behavior remapping
%% ============================================================

fprintf('classifier behavior -> canonical behavior mapping:\n');

for iBehv = 1:numel(classifierBehvs)

    canonicalInd = classBehvNumbers(iBehv);

    if canonicalInd == 0

        canonicalName = 'UNTRACKED / 0';

    else

        canonicalName = manBehvNames{canonicalInd};

    end

    fprintf('  %2d: %-20s -> %2d (%s)\n', ...
        iBehv, ...
        classifierBehvs{iBehv}, ...
        canonicalInd, ...
        canonicalName);

end

%% ============================================================
%% Build reverse map:
%% full downsampled-EMG index -> reduced behavior-label index
%% ============================================================

origDownsampEMGInd = origDownsampEMGInd(:);

if isempty(origDownsampEMGInd)
    error('origDownsampEMGInd is empty.');
end

maxFullIndex = max(origDownsampEMGInd);

map = nan(maxFullIndex, 1);

map(origDownsampEMGInd) = ...
    1:numel(origDownsampEMGInd);

%% ============================================================
%% Initialize outputs
%% ============================================================

nChannels = numel(validTransitionsCell);

regionLabelsPerTransition = cell(1, nChannels);
manualLabelsPerTransition = cell(1, nChannels);
classifierLabelsPerTransition = cell(1, nChannels);

%% ============================================================
%% Loop over each muscle channel
%% ============================================================

for ch = 1:nChannels

    validTransitions = validTransitionsCell{ch};

    regLabels = nan(size(validTransitions));
    manLabels = nan(size(validTransitions));
    classLabels = nan(size(validTransitions));

    for i = 1:numel(validTransitions)

        transitionIdx = validTransitions(i);

        %% ---------- map bounds ----------

        if transitionIdx < 1 || ...
                transitionIdx > numel(map)

            error( ...
                'mapTransitionsToCanonicalBehaviors:transitionIdxOutOfRange', ...
                ['transition index %d is out of range for map ' ...
                 '(1..%d)'], ...
                transitionIdx, ...
                numel(map));
        end

        %% ---------- convert to reduced behavior time base ----------

        umapIdx = map(transitionIdx);

        if isnan(umapIdx)

            % No corresponding reduced behavior-label sample.
            % Leave output labels as NaN.
            continue;

        end

        %% ========================================================
        %% Region label
        %% ========================================================
        % No separate old-style UMAP-region labeling exists for the
        % newer animals. Populate this compatibility field from the
        % classifier labels.

        regionVal = classifierLabels(umapIdx);

        if regionVal == 0

            regLabels(i) = 0;

        else

            if regionVal < 1 || ...
                    regionVal > numel(classBehvNumbers)

                error( ...
                    'mapTransitionsToCanonicalBehaviors:regionLabelOutOfRange', ...
                    ['region/classifier label index %d is out of ' ...
                     'range (1..%d)'], ...
                    regionVal, ...
                    numel(classBehvNumbers));
            end

            regLabels(i) = ...
                classBehvNumbers(regionVal);

        end

        %% ========================================================
        %% Manual label
        %% ========================================================
        % For these animals, populate the manual-label compatibility
        % field from classifierLabels.

        manualVal = classifierLabels(umapIdx);

        if manualVal == 0

            manLabels(i) = 0;

        else

            if manualVal < 1 || ...
                    manualVal > numel(manBehvNumbers)

                error( ...
                    'mapTransitionsToCanonicalBehaviors:manualLabelOutOfRange', ...
                    ['manual/classifier label index %d is out of ' ...
                     'range (1..%d)'], ...
                    manualVal, ...
                    numel(manBehvNumbers));
            end

            manLabels(i) = ...
                manBehvNumbers(manualVal);

        end

        %% ========================================================
        %% Classifier label
        %% ========================================================

        classVal = classifierLabels(umapIdx);

        if classVal == 0

            classLabels(i) = 0;

        else

            if classVal < 1 || ...
                    classVal > numel(classBehvNumbers)

                error( ...
                    'mapTransitionsToCanonicalBehaviors:classifierLabelOutOfRange', ...
                    ['classifier label index %d is out of range ' ...
                     '(1..%d)'], ...
                    classVal, ...
                    numel(classBehvNumbers));
            end

            classLabels(i) = ...
                classBehvNumbers(classVal);

        end
    end

    regionLabelsPerTransition{ch} = regLabels;
    manualLabelsPerTransition{ch} = manLabels;
    classifierLabelsPerTransition{ch} = classLabels;

    fprintf(['channel %d: %d transitions | ' ...
             '%d mapped | %d unmapped\n'], ...
        ch, ...
        numel(validTransitions), ...
        nnz(~isnan(classLabels)), ...
        nnz(isnan(classLabels)));
end

%% ============================================================
%% Save output variables
%% ============================================================

analyzedBehaviors = classifierBehvs;

outFile = fullfile( ...
    folderPath, ...
    'transitionCanonicalBehaviorLabels.mat');

save(outFile, ...
    'regionLabelsPerTransition', ...
    'manualLabelsPerTransition', ...
    'classifierLabelsPerTransition', ...
    'classifierBehvs', ...
    'manBehvNames', ...
    'analyzedBehaviors', ...
    'regionBehvAssignments', ...
    'mouseID', ...
    'baseSessionName', ...
    'probeRegion');

fprintf('\nsaved canonical transition labels:\n%s\n', outFile);

%% ============================================================
%% Plot histogram: region compatibility labels
%% ============================================================

figure( ...
    'Name', ...
    sprintf('%s Transitions by Region Label', mouseID), ...
    'Color', 'w');

tiledlayout(1, nChannels, ...
    'Padding', 'compact');

for ch = 1:nChannels

    nexttile;

    histogram( ...
        regionLabelsPerTransition{ch}, ...
        'BinMethod', 'integers', ...
        'FaceColor', 'b');

    xlabel('canonical region label (0-10)');
    ylabel('count');
    title(sprintf('muscle %d', ch));

end

%% ============================================================
%% Plot histogram: manual compatibility labels
%% ============================================================

figure( ...
    'Name', ...
    sprintf('%s Transitions by Manual Label', mouseID), ...
    'Color', 'w');

tiledlayout(1, nChannels, ...
    'Padding', 'compact');

for ch = 1:nChannels

    nexttile;

    histogram( ...
        manualLabelsPerTransition{ch}, ...
        'BinMethod', 'integers', ...
        'FaceColor', 'g');

    xlabel('manual canonical label (0-10)');
    ylabel('count');
    title(sprintf('muscle %d', ch));

end

%% ============================================================
%% Plot histogram: classifier canonical labels
%% ============================================================

figure( ...
    'Name', ...
    sprintf('%s Transitions by Classifier Label', mouseID), ...
    'Color', 'w');

tiledlayout(1, nChannels, ...
    'Padding', 'compact');

for ch = 1:nChannels

    nexttile;

    histogram( ...
        classifierLabelsPerTransition{ch}, ...
        'BinMethod', 'integers', ...
        'FaceColor', 'r');

    xlabel('classifier canonical label (0-10)');
    ylabel('count');
    title(sprintf('muscle %d', ch));

end

end
