function plotRawEMGDuringBehaviorAcrossAnimals(baseDirs, labelType, behaviorToPlot, maxEpochsPerAnimal)
% plots raw emg signal during all epochs of a chosen behavior across animals
% and marks detected emg transitions with red asterisks

% this is meant to check whether the emg during one behavior (for example,
% grooming = behavior 4) looks obviously different in one animal compared
% to the others

% what it does:
%   1. loads downsampEMG (1 kHz raw emg)
%   2. loads behavior labels aligned to the reduced time base from UMAP.mat
%   3. maps those reduced labels back onto the full 1 kHz emg timeline using
%      origDownsampEMGInd
%   4. finds contiguous epochs of the chosen behavior
%   5. plots raw emg snippets for each animal
%   6. marks detected transitions from validTransitionsCell as red asterisks

% inputs
%   baseDirs: cell array of session folders
%   labelType: currently only 'classifier' is supported cleanly
%   behaviorToPlot: canonical behavior number (default = 4 for grooming)
%   maxEpochsPerAnimal: maximum number of epochs to plot per animal (default = 20)

% j run:
% baseDirs = {
%     'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData', ...
%     'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData', ...
%     'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData', ...
%     'X:\David\ArenaRecordings\D043-020525-ArenaRecording\ProcessedData'
% };
% plotRawEMGDuringBehaviorAcrossAnimals(baseDirs, 'classifier', 4, 20)

if nargin < 2 || isempty(labelType)
    labelType = 'classifier';
end
if nargin < 3 || isempty(behaviorToPlot)
    behaviorToPlot = 4;
end
if nargin < 4 || isempty(maxEpochsPerAnimal)
    maxEpochsPerAnimal = 20;
end

labelType = lower(string(labelType));
if labelType ~= "classifier"
    error('this script is currently written for classifier labels.');
end

nAnimals = numel(baseDirs);
nCh = 4;

animalNames = cell(1, nAnimals);
allEpochStarts = cell(1, nAnimals);
allEpochEnds = cell(1, nAnimals);
allEpochDurations = cell(1, nAnimals);

for iDir = 1:nAnimals
    tok = regexp(baseDirs{iDir}, 'D\d+', 'match', 'once');
    if isempty(tok)
        animalNames{iDir} = sprintf('animal%d', iDir);
    else
        animalNames{iDir} = tok;
    end
end

%% ---------------- pass 1: build full-length behavior labels and epoch stats ----------------
for iDir = 1:nAnimals
    baseDir = baseDirs{iDir};

    E = load(fullfile(baseDir, 'EMG1ms.mat'), 'downsampEMG');
    U = load(fullfile(baseDir, 'UMAP.mat'), ...
        'origDownsampEMGInd', 'classifierLabels', 'classifierBehvs');

    nTime = size(E.downsampEMG, 2);

    labels1k = buildClassifierLabels1k(U, nTime);

    [epochStarts, epochEnds] = findBehaviorEpochs(labels1k, behaviorToPlot);

    allEpochStarts{iDir} = epochStarts;
    allEpochEnds{iDir} = epochEnds;
    allEpochDurations{iDir} = epochEnds - epochStarts + 1;
end

%% ---------------- summary printout ----------------
fprintf('\n=====================================================\n')
fprintf('behavior %d epoch summary (%s labels)\n', behaviorToPlot, labelType)
fprintf('=====================================================\n')

for iDir = 1:nAnimals
    durs = allEpochDurations{iDir};
    if isempty(durs)
        fprintf('%s: 0 epochs\n', animalNames{iDir});
    else
        fprintf('%s: %d epochs | median dur = %.1f ms | max dur = %d ms\n', ...
            animalNames{iDir}, numel(durs), median(durs), max(durs));
    end
end

%% ---------------- figure 1: first/longest raw epochs, stacked by animal ----------------
figure('Name', sprintf('Raw EMG during behavior %d', behaviorToPlot), 'Color', 'w');
tiledlayout(nAnimals, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

for iDir = 1:nAnimals
    baseDir = baseDirs{iDir};

    E = load(fullfile(baseDir, 'EMG1ms.mat'), 'downsampEMG');
    T = load(fullfile(baseDir, 'EMG_Neural_AllChannels.mat'), 'validTransitionsCell');
    emg = E.downsampEMG(1:nCh, :);
    validTransitionsCell = T.validTransitionsCell;

    epochStarts = allEpochStarts{iDir};
    epochEnds = allEpochEnds{iDir};
    epochDurations = allEpochDurations{iDir};

    nexttile
    hold on

    if isempty(epochStarts)
        title(sprintf('%s: no behavior-%d epochs', animalNames{iDir}, behaviorToPlot))
        xlabel('time within epoch (ms)')
        ylabel('emg')
        hold off
        continue
    end

    [~, sortIdx] = sort(epochDurations, 'descend');
    sortIdx = sortIdx(1:min(maxEpochsPerAnimal, numel(sortIdx)));

    offset = 0;
    offsetStep = 0;

    for k = 1:numel(sortIdx)
        idx = sortIdx(k);
        s = epochStarts(idx);
        e = epochEnds(idx);

        seg = emg(:, s:e);
        t = 0:(size(seg,2)-1);

        if k == 1
            offsetStep = max(abs(seg(:))) * 1.5;
            if offsetStep == 0 || isnan(offsetStep)
                offsetStep = 1;
            end
        end

        for ch = 1:nCh
            y = seg(ch,:) + offset + (ch-1)*offsetStep;
            plot(t, y, 'LineWidth', 1);

            % mark transitions from this channel that fall inside this epoch
            if ch <= numel(validTransitionsCell) && ~isempty(validTransitionsCell{ch})
                theseTransitions = validTransitionsCell{ch};
                inEpoch = theseTransitions(theseTransitions >= s & theseTransitions <= e);

                if ~isempty(inEpoch)
                    tRel = inEpoch - s;
                    yStars = emg(ch, inEpoch) + offset + (ch-1)*offsetStep;
                    plot(tRel, yStars, 'r*', 'MarkerSize', 5);
                end
            end
        end

        offset = offset + (nCh + 1) * offsetStep;
    end

    title(sprintf('%s: longest %d behavior-%d epochs', ...
        animalNames{iDir}, min(maxEpochsPerAnimal, numel(sortIdx)), behaviorToPlot))
    xlabel('time within epoch (ms)')
    ylabel('stacked raw emg')
    hold off
end

sgtitle(sprintf('behavior %d raw emg with transition markers', behaviorToPlot))

%% ---------------- figure 2: overlay mean absolute emg envelope during behavior ----------------
figure('Name', sprintf('Mean abs EMG during behavior %d', behaviorToPlot), 'Color', 'w');
tiledlayout(1, nCh, 'Padding', 'compact', 'TileSpacing', 'compact');

allDurVec = vertcat(allEpochDurations{:});
if isempty(allDurVec)
    targetLen = 1;
else
    targetLen = max(1, round(median(allDurVec)));
end

for ch = 1:nCh
    nexttile
    hold on

    for iDir = 1:nAnimals
        baseDir = baseDirs{iDir};

        E = load(fullfile(baseDir, 'EMG1ms.mat'), 'downsampEMG');
        emg = E.downsampEMG(ch, :);

        epochStarts = allEpochStarts{iDir};
        epochEnds = allEpochEnds{iDir};

        if isempty(epochStarts)
            continue
        end

        X = nan(numel(epochStarts), targetLen);

        for k = 1:numel(epochStarts)
            s = epochStarts(k);
            e = epochEnds(k);

            seg = abs(emg(s:e));
            L = min(numel(seg), targetLen);
            X(k,1:L) = seg(1:L);
        end

        mu = mean(X, 1, 'omitnan');
        plot(0:(targetLen-1), mu, 'LineWidth', 1.5, 'DisplayName', animalNames{iDir});
    end

    title(sprintf('channel %d', ch))
    xlabel('time from epoch start (ms)')
    ylabel('mean |emg|')
    legend('Location', 'best')
    hold off
end

%% ---------------- figure 3: zoomed individual epochs for each animal/channel ----------------
figure('Name', sprintf('Behavior %d raw EMG snippets by animal', behaviorToPlot), 'Color', 'w');
tiledlayout(nAnimals, nCh, 'Padding', 'compact', 'TileSpacing', 'compact');

for iDir = 1:nAnimals
    baseDir = baseDirs{iDir};

    E = load(fullfile(baseDir, 'EMG1ms.mat'), 'downsampEMG');
    T = load(fullfile(baseDir, 'EMG_Neural_AllChannels.mat'), 'validTransitionsCell');
    emg = E.downsampEMG(1:nCh, :);
    validTransitionsCell = T.validTransitionsCell;

    epochStarts = allEpochStarts{iDir};
    epochEnds = allEpochEnds{iDir};
    epochDurations = allEpochDurations{iDir};

    if isempty(epochStarts)
        for ch = 1:nCh
            nexttile
            title(sprintf('%s ch %d (n=0)', animalNames{iDir}, ch))
            xlabel('time (ms)')
            ylabel('emg')
        end
        continue
    end

    [~, sortIdx] = sort(epochDurations, 'descend');
    sortIdx = sortIdx(1:min(maxEpochsPerAnimal, numel(sortIdx)));

    for ch = 1:nCh
        nexttile
        hold on

        for k = 1:numel(sortIdx)
            idx = sortIdx(k);
            s = epochStarts(idx);
            e = epochEnds(idx);

            seg = emg(ch, s:e);
            t = 0:(numel(seg)-1);
            plot(t, seg, 'Color', [0.7 0.7 0.7]);

            % mark transitions in this epoch for this channel
            if ch <= numel(validTransitionsCell) && ~isempty(validTransitionsCell{ch})
                theseTransitions = validTransitionsCell{ch};
                inEpoch = theseTransitions(theseTransitions >= s & theseTransitions <= e);

                if ~isempty(inEpoch)
                    tRel = inEpoch - s;
                    yStars = emg(ch, inEpoch);
                    plot(tRel, yStars, 'r*', 'MarkerSize', 5);
                end
            end
        end

        title(sprintf('%s ch %d (n=%d)', animalNames{iDir}, ch, numel(sortIdx)))
        xlabel('time within epoch (ms)')
        ylabel('raw emg')
        hold off
    end
end

sgtitle(sprintf('behavior %d raw emg snippets with transition markers', behaviorToPlot))

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

labelVecReduced = nan(size(classifierLabels));
for i = 1:numel(classifierLabels)
    v = classifierLabels(i);

    if isnan(v)
        continue
    elseif v == 0
        labelVecReduced(i) = 0;
    else
        if v >= 1 && v <= numel(classBehvNumbers)
            labelVecReduced(i) = classBehvNumbers(v);
        else
            labelVecReduced(i) = nan;
        end
    end
end

n = min(numel(origInd), numel(labelVecReduced));
origInd = origInd(1:n);
labelVecReduced = labelVecReduced(1:n);

ok = origInd >= 1 & origInd <= nTime & ~isnan(labelVecReduced);
labels1k(origInd(ok)) = labelVecReduced(ok);

end

%% ============================================================
function [epochStarts, epochEnds] = findBehaviorEpochs(labels1k, behaviorToPlot)
mask = labels1k == behaviorToPlot;
mask(isnan(mask)) = false;

d = diff([false, mask, false]);
epochStarts = find(d == 1);
epochEnds = find(d == -1) - 1;

end
