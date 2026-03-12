function plotBehavior4EMGWithRecomputedThresholdCrossings(baseDirs, threshold, baselineDur, minSeparation, channelToPlot, maxEpochsPerAnimal)
% plots raw emg during all behavior-4 epochs and recomputes threshold crossings
% directly from the raw signal, marking those crossings with red asterisks

% this version supports plotting multiple channels, e.g. channelToPlot = 1:4

% inputs
%   baseDirs : cell array of processeddata folders
%   threshold : emg threshold for defining activation crossings
%   baselineDur : minimum quiet baseline before a transition (ms)
%   minSeparation : minimum separation between crossings (ms)
%   channelToPlot : scalar or vector of emg channels to plot (default = 1)
%   maxEpochsPerAnimal: maximum number of behavior-4 epochs to plot (default = 20)

% transition logic matches the earlier extraction logic:
%   1. signal > threshold
%   2. take rising edges only
%   3. enforce min separation
%   4. require quiet baseline for baselineDur before the crossing

% j run:
% baseDirs = {
%     'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData', ...
%     'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData', ...
%     'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData', ...
%     'X:\David\ArenaRecordings\D043-020525-ArenaRecording\ProcessedData'
% };
% plotBehavior4EMGWithRecomputedThresholdCrossings(baseDirs, 100, 500, 500, 1:4, 20)

if nargin < 2 || isempty(threshold)
    threshold = 100;
end
if nargin < 3 || isempty(baselineDur)
    baselineDur = 500;
end
if nargin < 4 || isempty(minSeparation)
    minSeparation = 500;
end
if nargin < 5 || isempty(channelToPlot)
    channelToPlot = 1;
end
if nargin < 6 || isempty(maxEpochsPerAnimal)
    maxEpochsPerAnimal = 20;
end

behaviorToPlot = 4;
channelsToPlot = channelToPlot(:)';
nChPlot = numel(channelsToPlot);
nAnimals = numel(baseDirs);

animalNames = cell(1, nAnimals);
for iDir = 1:nAnimals
    tok = regexp(baseDirs{iDir}, 'D\d+', 'match', 'once');
    if isempty(tok)
        animalNames{iDir} = sprintf('animal%d', iDir);
    else
        animalNames{iDir} = tok;
    end
end

allEpochStarts = cell(1, nAnimals);
allEpochEnds = cell(1, nAnimals);
allEpochDurations = cell(1, nAnimals);

% store transitions per animal per plotted channel
allDetectedTransitions = cell(nAnimals, nChPlot);
allTransitionsInBehavior4 = cell(nAnimals, nChPlot);

%% ---------------- pass 1: build behavior-4 epochs and recompute crossings ----------------
for iDir = 1:nAnimals
    baseDir = baseDirs{iDir};

    E = load(fullfile(baseDir, 'EMG1ms.mat'), 'downsampEMG');
    downsampEMG = E.downsampEMG;

    allSignal = downsampEMG(1:4, :);
    nTime = size(downsampEMG, 2);

    U = load(fullfile(baseDir, 'UMAP.mat'), ...
        'origDownsampEMGInd', 'classifierLabels', 'classifierBehvs');

    labels1k = buildClassifierLabels1k(U, nTime);

    [epochStarts, epochEnds] = findBehaviorEpochs(labels1k, behaviorToPlot);
    epochDurations = epochEnds - epochStarts + 1;

    allEpochStarts{iDir} = epochStarts;
    allEpochEnds{iDir} = epochEnds;
    allEpochDurations{iDir} = epochDurations;

    for iCh = 1:nChPlot
        ch = channelsToPlot(iCh);
        signal = downsampEMG(ch, :);

        above = signal > threshold;
        transitionPoints = find(diff([0, above]) == 1);

        if numel(transitionPoints) > 1
            diffs = diff(transitionPoints);
            keepers = [true, diffs > minSeparation];
            transitionPoints = transitionPoints(keepers);
        end

        validTransitions = [];
        for i = 1:numel(transitionPoints)
            idx = transitionPoints(i);

            if idx > baselineDur && all(allSignal(:, idx-baselineDur:idx-1) < threshold, 'all')
                validTransitions(end+1) = idx; %#ok<AGROW>
            end
        end

        allDetectedTransitions{iDir, iCh} = validTransitions;

        transInBeh4 = [];
        for k = 1:numel(epochStarts)
            s = epochStarts(k);
            e = epochEnds(k);
            transInBeh4 = [transInBeh4; validTransitions(validTransitions >= s & validTransitions <= e)']; %#ok<AGROW>
        end
        allTransitionsInBehavior4{iDir, iCh} = unique(transInBeh4(:))';
    end
end

%% ---------------- summary printout ----------------
fprintf('\n=====================================================\n')
fprintf('behavior 4 epoch summary with recomputed threshold crossings\n')
fprintf('=====================================================\n')
fprintf('threshold = %g | baselineDur = %d ms | minSeparation = %d ms\n\n', ...
    threshold, baselineDur, minSeparation)

for iDir = 1:nAnimals
    durs = allEpochDurations{iDir};

    if isempty(durs)
        fprintf('%s: 0 behavior-4 epochs\n', animalNames{iDir});
    else
        fprintf('%s: %d behavior-4 epochs | median dur = %.1f ms | max dur = %d ms\n', ...
            animalNames{iDir}, numel(durs), median(durs), max(durs));
    end

    for iCh = 1:nChPlot
        ch = channelsToPlot(iCh);
        fprintf('    ch %d | total crossings = %d | behavior-4 crossings = %d\n', ...
            ch, numel(allDetectedTransitions{iDir, iCh}), numel(allTransitionsInBehavior4{iDir, iCh}));
    end
end

%% ---------------- figure 1: longest behavior-4 epochs with red asterisks ----------------
figure('Name', 'behavior 4 raw emg with recomputed threshold crossings', 'Color', 'w');
tiledlayout(nAnimals, nChPlot, 'Padding', 'compact', 'TileSpacing', 'compact');

for iDir = 1:nAnimals
    baseDir = baseDirs{iDir};

    E = load(fullfile(baseDir, 'EMG1ms.mat'), 'downsampEMG');
    downsampEMG = E.downsampEMG;

    epochStarts = allEpochStarts{iDir};
    epochEnds = allEpochEnds{iDir};
    epochDurations = allEpochDurations{iDir};

    for iCh = 1:nChPlot
        ch = channelsToPlot(iCh);
        signal = downsampEMG(ch, :);
        validTransitions = allDetectedTransitions{iDir, iCh};

        nexttile
        hold on

        if isempty(epochStarts)
            title(sprintf('%s ch %d: no behavior-4 epochs', animalNames{iDir}, ch))
            xlabel('time within epoch (ms)')
            ylabel('raw emg')
            hold off
            continue
        end

        [~, sortIdx] = sort(epochDurations, 'descend');
        sortIdx = sortIdx(1:min(maxEpochsPerAnimal, numel(sortIdx)));

        for k = 1:numel(sortIdx)
            idx = sortIdx(k);
            s = epochStarts(idx);
            e = epochEnds(idx);

            seg = signal(s:e);
            t = 0:(numel(seg)-1);

            plot(t, seg, 'Color', [0.7 0.7 0.7]);

            inEpoch = validTransitions(validTransitions >= s & validTransitions <= e);
            if ~isempty(inEpoch)
                tRel = inEpoch - s;
                yStar = signal(inEpoch);
                plot(tRel, yStar, 'r*', 'MarkerSize', 5);
            end
        end

        title(sprintf('%s ch %d (n=%d)', ...
            animalNames{iDir}, ch, min(maxEpochsPerAnimal, numel(sortIdx))))
        xlabel('time within epoch (ms)')
        ylabel('raw emg')
        hold off
    end
end

sgtitle(sprintf('behavior 4 raw emg with recomputed threshold crossings (threshold = %g)', threshold))

%% ---------------- figure 2: stacked epochs with red asterisks ----------------
figure('Name', 'stacked behavior 4 raw emg with recomputed crossings', 'Color', 'w');
tiledlayout(nAnimals, nChPlot, 'Padding', 'compact', 'TileSpacing', 'compact');

for iDir = 1:nAnimals
    baseDir = baseDirs{iDir};

    E = load(fullfile(baseDir, 'EMG1ms.mat'), 'downsampEMG');
    downsampEMG = E.downsampEMG;

    epochStarts = allEpochStarts{iDir};
    epochEnds = allEpochEnds{iDir};
    epochDurations = allEpochDurations{iDir};

    for iCh = 1:nChPlot
        ch = channelsToPlot(iCh);
        signal = downsampEMG(ch, :);
        validTransitions = allDetectedTransitions{iDir, iCh};

        nexttile
        hold on

        if isempty(epochStarts)
            title(sprintf('%s ch %d: no behavior-4 epochs', animalNames{iDir}, ch))
            xlabel('time within epoch (ms)')
            ylabel('stacked raw emg')
            hold off
            continue
        end

        [~, sortIdx] = sort(epochDurations, 'descend');
        sortIdx = sortIdx(1:min(maxEpochsPerAnimal, numel(sortIdx)));

        offset = 0;
        offsetStep = max(abs(signal(:))) * 1.2;
        if offsetStep == 0 || isnan(offsetStep)
            offsetStep = 1;
        end

        for k = 1:numel(sortIdx)
            idx = sortIdx(k);
            s = epochStarts(idx);
            e = epochEnds(idx);

            seg = signal(s:e);
            t = 0:(numel(seg)-1);

            plot(t, seg + offset, 'k')

            inEpoch = validTransitions(validTransitions >= s & validTransitions <= e);
            if ~isempty(inEpoch)
                tRel = inEpoch - s;
                yStar = signal(inEpoch) + offset;
                plot(tRel, yStar, 'r*', 'MarkerSize', 5);
            end

            offset = offset + offsetStep;
        end

        title(sprintf('%s ch %d (n=%d)', ...
            animalNames{iDir}, ch, min(maxEpochsPerAnimal, numel(sortIdx))))
        xlabel('time within epoch (ms)')
        ylabel('stacked raw emg')
        hold off
    end
end

sgtitle('behavior 4 raw emg with recomputed threshold crossings')

end

%% ============================================================
function labels1k = buildClassifierLabels1k(U, nTime)
origInd = U.origDownsampEMGInd(:);
classifierLabels = U.classifierLabels(:);
classifierBehvs = U.classifierBehvs;

labels1k = nan(1, nTime);

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

nClassBehv = numel(classifierBehvs);
classBehvNumbers = zeros(1, nClassBehv);

for iBehv = 1:nClassBehv
    thisName = classifierBehvs{iBehv};
    cleanName = lower(strrep(strrep(thisName, ' ', ''), '_', ''));

    if isKey(canonicalLookup, cleanName)
        classBehvNumbers(iBehv) = canonicalLookup(cleanName);
    else
        classBehvNumbers(iBehv) = 0;
    end
end

n = min(numel(origInd), numel(classifierLabels));
origInd = origInd(1:n);
classifierLabels = classifierLabels(1:n);

for i = 1:n
    idx = classifierLabels(i);

    if isnan(idx) || idx < 0
        continue
    elseif idx == 0
        val = 0;
    else
        if idx <= numel(classBehvNumbers)
            val = classBehvNumbers(idx);
        else
            continue
        end
    end

    emgIdx = origInd(i);
    if emgIdx >= 1 && emgIdx <= nTime
        labels1k(emgIdx) = val;
    end
end

end

%% ============================================================
function [epochStarts, epochEnds] = findBehaviorEpochs(labels1k, behaviorToPlot)
mask = labels1k == behaviorToPlot;
mask(isnan(mask)) = false;

d = diff([false, mask, false]);
epochStarts = find(d == 1);
epochEnds = find(d == -1) - 1;

end
