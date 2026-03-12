function plotGroomingEMGAcrossAnimals(baseDirs, labelType, behaviorToPlot)
% plots emg windows for one behavior across animals
% intended for checking whether grooming-labeled emg looks odd in one animal

% this uses transition-level labels already matched to emg/neural windows and plots emg around those labeled transitions for each muscle channel

% inputs
%   baseDirs : cell array of ProcessedData/session folders
%   labelType : 'classifier', 'manual', or 'umap' (default = 'classifier')
%   behaviorToPlot : canonical behavior index to plot (default = 4 for grooming)
%
% j run:
% baseDirs = {
%     'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData', ...
%     'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData', ...
%     'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData', ...
%     'X:\David\ArenaRecordings\D043-020525-ArenaRecording\ProcessedData'
% };
% plotGroomingEMGAcrossAnimals(baseDirs, 'classifier', 4)

if nargin < 2 || isempty(labelType)
    labelType = 'classifier';
end
if nargin < 3 || isempty(behaviorToPlot)
    behaviorToPlot = 4;
end

labelType = lower(string(labelType));
nAnimals = numel(baseDirs);
nCh = 4;

% animal names for titles
animalNames = cell(1, nAnimals);
for iDir = 1:nAnimals
    thisDir = baseDirs{iDir};
    tok = regexp(thisDir, 'D\d+', 'match', 'once');
    if isempty(tok)
        animalNames{iDir} = sprintf('animal %d', iDir);
    else
        animalNames{iDir} = tok;
    end
end

% store counts for a summary printout
trialCountMat = zeros(nAnimals, nCh);

% ---------------- figure 1: mean ± sem emg traces ----------------
figure('Name', sprintf('behavior %d emg mean traces across animals', behaviorToPlot), 'Color', 'w');
tiledlayout(nAnimals, nCh, 'Padding', 'compact', 'TileSpacing', 'compact');

for iDir = 1:nAnimals
    baseDir = baseDirs{iDir};

    % load emg windows
    S = load(fullfile(baseDir, 'EMG_Neural_AllChannels.mat'), ...
        'emgWindowsCell', 'tAxis');

    % load labels
    L = load(fullfile(baseDir, 'transitionCanonicalBehaviorLabels.mat'), ...
        'regionLabelsPerTransition', 'manualLabelsPerTransition', 'classifierLabelsPerTransition');

    switch labelType
        case "umap"
            labelsPerTransition = L.regionLabelsPerTransition;
        case "manual"
            labelsPerTransition = L.manualLabelsPerTransition;
        case "classifier"
            labelsPerTransition = L.classifierLabelsPerTransition;
        otherwise
            error('unknown labelType: %s', labelType);
    end

    tAxis = S.tAxis(:)';

    for ch = 1:nCh
        nexttile

        if ch > numel(S.emgWindowsCell) || isempty(S.emgWindowsCell{ch}) || ...
                ch > numel(labelsPerTransition) || isempty(labelsPerTransition{ch})
            title(sprintf('%s ch %d (empty)', animalNames{iDir}, ch))
            xlabel('time (ms)')
            ylabel('emg')
            continue
        end

        emgWin = S.emgWindowsCell{ch};      % events x muscles x time
        labels = labelsPerTransition{ch}(:);

        nEvt = min(size(emgWin,1), numel(labels));
        emgWin = emgWin(1:nEvt,:,:);
        labels = labels(1:nEvt);

        keepIdx = find(labels == behaviorToPlot);
        trialCountMat(iDir, ch) = numel(keepIdx);

        if isempty(keepIdx)
            title(sprintf('%s ch %d (n=0)', animalNames{iDir}, ch))
            xlabel('time (ms)')
            ylabel('emg')
            xlim([tAxis(1) tAxis(end)])
            continue
        end

        % use the same muscle channel as the transition channel for the main comparison
        thisEMG = squeeze(emgWin(keepIdx, ch, :));   % trials x time

        if isvector(thisEMG)
            thisEMG = reshape(thisEMG, 1, []);
        end

        mu = mean(thisEMG, 1, 'omitnan');
        sem = std(thisEMG, 0, 1, 'omitnan') ./ sqrt(size(thisEMG,1));

        hold on
        fill([tAxis fliplr(tAxis)], [mu-sem fliplr(mu+sem)], ...
            [0.8 0.85 1], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        plot(tAxis, mu, 'b', 'LineWidth', 1.5);
        xline(0, '--k');
        hold off

        title(sprintf('%s ch %d (n=%d)', animalNames{iDir}, ch, numel(keepIdx)))
        xlabel('time (ms)')
        ylabel('emg')
        xlim([tAxis(1) tAxis(end)])
    end
end

sgtitle(sprintf('behavior %d emg mean ± sem by animal and muscle', behaviorToPlot))

% ---------------- figure 2: single-trial traces ----------------
figure('Name', sprintf('behavior %d emg single trials across animals', behaviorToPlot), 'Color', 'w');
tiledlayout(nAnimals, nCh, 'Padding', 'compact', 'TileSpacing', 'compact');

for iDir = 1:nAnimals
    baseDir = baseDirs{iDir};

    S = load(fullfile(baseDir, 'EMG_Neural_AllChannels.mat'), ...
        'emgWindowsCell', 'tAxis');

    L = load(fullfile(baseDir, 'transitionCanonicalBehaviorLabels.mat'), ...
        'regionLabelsPerTransition', 'manualLabelsPerTransition', 'classifierLabelsPerTransition');

    switch labelType
        case "umap"
            labelsPerTransition = L.regionLabelsPerTransition;
        case "manual"
            labelsPerTransition = L.manualLabelsPerTransition;
        case "classifier"
            labelsPerTransition = L.classifierLabelsPerTransition;
    end

    tAxis = S.tAxis(:)';

    for ch = 1:nCh
        nexttile

        if ch > numel(S.emgWindowsCell) || isempty(S.emgWindowsCell{ch}) || ...
                ch > numel(labelsPerTransition) || isempty(labelsPerTransition{ch})
            title(sprintf('%s ch %d (empty)', animalNames{iDir}, ch))
            xlabel('time (ms)')
            ylabel('emg')
            continue
        end

        emgWin = S.emgWindowsCell{ch};
        labels = labelsPerTransition{ch}(:);

        nEvt = min(size(emgWin,1), numel(labels));
        emgWin = emgWin(1:nEvt,:,:);
        labels = labels(1:nEvt);

        keepIdx = find(labels == behaviorToPlot);

        if isempty(keepIdx)
            title(sprintf('%s ch %d (n=0)', animalNames{iDir}, ch))
            xlabel('time (ms)')
            ylabel('emg')
            xlim([tAxis(1) tAxis(end)])
            continue
        end

        thisEMG = squeeze(emgWin(keepIdx, ch, :));   % trials x time

        if isvector(thisEMG)
            thisEMG = reshape(thisEMG, 1, []);
        end

        hold on
        for tr = 1:size(thisEMG,1)
            plot(tAxis, thisEMG(tr,:), 'Color', [0.6 0.6 0.6]);
        end
        plot(tAxis, mean(thisEMG, 1, 'omitnan'), 'b', 'LineWidth', 1.5);
        xline(0, '--k');
        hold off

        title(sprintf('%s ch %d (n=%d)', animalNames{iDir}, ch, numel(keepIdx)))
        xlabel('time (ms)')
        ylabel('emg')
        xlim([tAxis(1) tAxis(end)])
    end
end

sgtitle(sprintf('behavior %d emg single trials by animal and muscle', behaviorToPlot))

% ---------------- print counts ----------------
fprintf('\n=============================================\n')
fprintf('behavior %d trial counts per animal and muscle\n', behaviorToPlot)
fprintf('=============================================\n')

for iDir = 1:nAnimals
    fprintf('%s: ', animalNames{iDir})
    fprintf('ch1=%d  ch2=%d  ch3=%d  ch4=%d\n', ...
        trialCountMat(iDir,1), trialCountMat(iDir,2), ...
        trialCountMat(iDir,3), trialCountMat(iDir,4))
end

fprintf('\n')
end
