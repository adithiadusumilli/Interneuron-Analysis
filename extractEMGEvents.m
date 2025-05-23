function extractEMGEvents(emgFile, threshold, baselineDur, minSeparation)
% function to identify threshold-crossing activation events in multi-channel emg and extract signal windows around those points

% inputs:
%   emgfile       — path to .mat file containing downsampemg (8 × timepoints)
%   threshold     — activation threshold value (e.g., 100)
%   baselinedur   — minimum quiet period (in ms) before activation (default: 200)
%   minseparation — minimum ms between separate events (default: 500)

    if nargin < 3, baselineDur = 200; end
    if nargin < 4, minSeparation = 500; end

    load(emgFile, 'downsampEMG');  % must contain 8 × n matrix

    [nChannels, nPoints] = size(downsampEMG);
    fs = 1000;  % assume 1 khz sampling rate (1 sample = 1 ms)

    preSamples = 100;   % number of samples before threshold crossing to include in window
    postSamples = 100;  % number of samples after threshold crossing to include in window

    for ch = 1:nChannels
        signal = downsampEMG(ch, :);  % get emg signal for current channel

        % apply threshold to find where signal is above threshold
        above = signal > threshold;

        % find indices where signal rises from below to above threshold
        transitionPoints = find(diff([0, above]) == 1);

        % remove transitions that are too close together (< minseparation)
        if numel(transitionPoints) > 1
            diffs = diff(transitionPoints);
            keepers = [true, diffs > minSeparation]; % bool for filtering the close proximity threshold crossings
            transitionPoints = transitionPoints(keepers);
        end

        % filter out threshold crossings that do not have quiet baseline before
        validTransitions = [];
        for i = 1:numel(transitionPoints)
            idx = transitionPoints(i);
            if idx > baselineDur && all(signal(idx-baselineDur:idx-1) < threshold)
                validTransitions(end+1) = idx;  % saving clean transition point
            end
        end

        % extract 201-sample window (100 before, 100 after, plus crossing) around each valid transition
        windows = [];
        for i = 1:numel(validTransitions)
            idx = validTransitions(i);
            if idx > preSamples && idx + postSamples <= nPoints
                snippet = signal(idx-preSamples:idx+postSamples);  % extract emg window
                windows(end+1, :) = snippet;  % add to array
            end
        end

        % plot emg signal and mark detected threshold crossings with red stars
        figure('Name', sprintf('emg channel %d', ch));
        plot(signal);
        hold on;
        plot(validTransitions, signal(validTransitions), 'r*');
        title(sprintf('channel %d emg with detected transitions', ch));
        xlabel('time (ms)');
        ylabel('emg amplitude');
        hold off;

        % save results for this channel
        saveFile = sprintf('EMG_Channel%d_Events.mat', ch);
        save(saveFile, 'validTransitions', 'windows');
        fprintf('channel %d: saved %d events to %s\n', ch, size(windows, 1), saveFile);
    end
end