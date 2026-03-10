function mapTransitionsToCanonicalBehaviors(folderPath)
% map emg transition windows to canonical behavior labels from umap, manual annotation, and classifier-based labels
% no actual umap or manual labels, they are just copies of classifier labeling (new animals dont have manual and UMAP format is diff)

% this function:
%   - loads emg transition indices for each muscle channel
%   - maps each transition index into the reduced time axis used for umap / manual labels
%   - assigns behavior labels using a consistent canonical scheme across animals:
%         umap regions:          1–7                       → regionLabelsPerTransition
%         manual behaviors:      0–10 (0 = unlabeled)      → manualLabelsPerTransition
%         classifier behaviors:  0–10 (0 = unlabeled)      → classifierLabelsPerTransition
%   - manual and classifier labels are remapped so that:
%         1..10 always correspond to the same ordered list in manBehvNames
%         any behavior not in manBehvNames (e.g. nosepoke) is mapped to 0
%   - saves outputs and plots per-muscle histograms for each label type

% load transition indices (in downsampemg time base)
load(fullfile(folderPath, 'EMG_Neural_AllChannels.mat'), 'validTransitionsCell');

% load behavior variables from umap
load(fullfile(folderPath, 'UMAP.mat'), ...
    'origDownsampEMGInd', ...          % mapping from reduced umap index to full downsampemg index
    'classifierLabels', ...            % classifier label indices (0 = unlabeled)
    'classifierBehvs', ...             % behavior names for classifier labels (per animal)
    'regionBehvAssignments');          % kept for saving to preserve pipeline compatibility

% canonical behavior name list and ordering (shared across animals)
% indices 1..10 in the canonical space will always correspond to these names
manBehvNames = {'climbdown','climbup','eating','grooming', ...
                'jumpdown','jumping','rearing','still','walkflat','walkgrid'};

% lookup table for alternate behavior names across animals / label sources
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

%% build manual-label remap: per-animal manual index -> canonical index (0..10)

nManualBehv = numel(classifierBehvs);
manBehvNumbers = zeros(1, nManualBehv);  % 0 = behavior not in canonical list (e.g. nosepoke)

for iBehv = 1:nManualBehv
    thisName = classifierBehvs{iBehv};
    cleanName = lower(strrep(strrep(thisName, ' ', ''), '_', ''));

    if isKey(canonicalLookup, cleanName)
        manBehvNumbers(iBehv) = canonicalLookup(cleanName);
    else
        % this behavior name is not in the canonical list
        % thus map to 0 and treat as unlabeled in canonical manual-label space
        manBehvNumbers(iBehv) = 0;
    end
end

%% build classifier-label remap: per-animal classifier index -> canonical index (0..10)

nClassBehv = numel(classifierBehvs);
classBehvNumbers = zeros(1, nClassBehv);  % 0 = behavior not in canonical list

for iBehv = 1:nClassBehv
    thisName = classifierBehvs{iBehv};
    cleanName = lower(strrep(strrep(thisName, ' ', ''), '_', ''));

    if isKey(canonicalLookup, cleanName)
        classBehvNumbers(iBehv) = canonicalLookup(cleanName);
    else
        % classifier behavior not in canonical list → map to 0
        classBehvNumbers(iBehv) = 0;
    end
end

%% build reverse map from full (~5m) downsampemg index to reduced umap index (~4m)

% map(k) gives the umap index corresponding to downsampemg index k
map = nan(max(origDownsampEMGInd), 1);
map(origDownsampEMGInd) = 1:numel(origDownsampEMGInd);

% outputs: 1 cell per muscle (1x4), each storing label per transition
regionLabelsPerTransition = cell(1, 4);
manualLabelsPerTransition = cell(1, 4);
classifierLabelsPerTransition = cell(1, 4);

%% loop over each muscle channel (1–4)

for ch = 1:4
    validTransitions = validTransitionsCell{ch};  % transition indices for this muscle

    regLabels = nan(size(validTransitions));    % umap region labels 1–7
    manLabels = nan(size(validTransitions));    % manual canonical labels 0–10
    classLabels = nan(size(validTransitions));  % classifier canonical labels 0–10

    for i = 1:length(validTransitions)
        transitionIdx = validTransitions(i);

        % sanity check for map bounds
        if transitionIdx < 1 || transitionIdx > numel(map)
            error('mapTransitionsToCanonicalBehaviors:transitionIdxOutOfRange', ...
                'transition index %d out of range for map (1..%d)', ...
                transitionIdx, numel(map));
        end

        % get the index into the reduced umap / behavior-label time base
        umapIdx = map(transitionIdx);

        if isnan(umapIdx)
            % this transition does not have a corresponding umap/manual index
            % leave all labels as nan and continue
            continue;
        end

        %% umap region mapping → use classifier labels to preserve pipeline
        % region labels no longer exist separately, so keep the field
        % but populate it from classifier labels

        regionVal = classifierLabels(umapIdx);  % 0 = unlabeled, >0 indexes classifierBehvs

        if regionVal == 0
            regLabels(i) = 0;
        else
            if regionVal < 1 || regionVal > numel(classBehvNumbers)
                error('mapTransitionsToCanonicalBehaviors:regionLabelOutOfRange', ...
                    'region label index %d out of range (1..%d)', ...
                    regionVal, numel(classBehvNumbers));
            end
            regLabels(i) = classBehvNumbers(regionVal);
        end

        %% manual label mapping → use classifier labels to preserve pipeline
        % manual labels no longer exist separately, so keep the manual field
        % but populate it using classifier labels to preserve downstream code

        manualVal = classifierLabels(umapIdx);  % 0 = unlabeled, >0 indexes classifierBehvs

        if manualVal == 0
            % unlabeled in original annotations → keep as 0
            manLabels(i) = 0;
        else
            if manualVal < 1 || manualVal > numel(manBehvNumbers)
                error('mapTransitionsToCanonicalBehaviors:manualLabelOutOfRange', ...
                    'manual label index %d out of range (1..%d)', ...
                    manualVal, numel(manBehvNumbers));
            end
            % remap to canonical index; can be 0 if this behavior is not tracked
            manLabels(i) = manBehvNumbers(manualVal);
        end

        %% classifier label mapping → canonical 0..10

        classVal = classifierLabels(umapIdx);  % 0 = unlabeled, >0 indexes classifierBehvs

        if classVal == 0
            % unlabeled by classifier → keep as 0
            classLabels(i) = 0;
        else
            if classVal < 1 || classVal > numel(classBehvNumbers)
                error('mapTransitionsToCanonicalBehaviors:classifierLabelOutOfRange', ...
                    'classifier label index %d out of range (1..%d)', ...
                    classVal, numel(classBehvNumbers));
            end
            % remap to canonical index; can be 0 if this classifier behavior is not tracked
            classLabels(i) = classBehvNumbers(classVal);
        end
    end

    regionLabelsPerTransition{ch} = regLabels;
    manualLabelsPerTransition{ch} = manLabels;
    classifierLabelsPerTransition{ch} = classLabels;
end

%% save output variables (including canonical behavior name list)

analyzedBehaviors = classifierBehvs;

save(fullfile(folderPath, 'transitionCanonicalBehaviorLabels.mat'), ...
    'regionLabelsPerTransition', 'manualLabelsPerTransition', ...
    'classifierLabelsPerTransition', ...
    'classifierBehvs', 'manBehvNames', 'analyzedBehaviors', ...
    'regionBehvAssignments');

%% plot histogram of transitions by umap region for each muscle

figure('Name', 'Transitions by UMAP Region (Per Muscle)', 'Color', 'w');
tiledlayout(1, 4, 'Padding', 'compact');
for ch = 1:4
    nexttile;
    histogram(regionLabelsPerTransition{ch}, 'BinMethod', 'integers', 'FaceColor', 'b');
    xlabel('umap region (1–7)');
    ylabel('count');
    title(sprintf('muscle %d', ch));
end

%% plot histogram of transitions by manual canonical label for each muscle

figure('Name', 'Transitions by Manual Label (Per Muscle)', 'Color', 'w');
tiledlayout(1, 4, 'Padding', 'compact');
for ch = 1:4
    nexttile;
    histogram(manualLabelsPerTransition{ch}, 'BinMethod', 'integers', 'FaceColor', 'g');
    xlabel('manual canonical label (0–10)');
    ylabel('count');
    title(sprintf('muscle %d', ch));
end

%% plot histogram of transitions by classifier canonical label for each muscle

figure('Name', 'Transitions by Classifier Label (Per Muscle)', 'Color', 'w');
tiledlayout(1, 4, 'Padding', 'compact');
for ch = 1:4
    nexttile;
    histogram(classifierLabelsPerTransition{ch}, 'BinMethod', 'integers', 'FaceColor', 'r');
    xlabel('classifier canonical label (0–10)');
    ylabel('count');
    title(sprintf('muscle %d', ch));
end
end
