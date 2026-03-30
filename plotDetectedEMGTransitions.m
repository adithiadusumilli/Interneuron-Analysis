function plotDetectedEMGTransitions(emgNeuralFile, emgChannel, t0, t1)
% plots one emg channel over a selected time window and overlays the
% already-detected valid transition points saved in EMG_Neural_AllChannels.mat

% inputs:
%   emgNeuralFile: full file path to EMG_Neural_AllChannels.mat
%   emgChannel: emg channel to plot
%   t0: start time in seconds
%   t1: end time in seconds

% this function:
%   1. loads emg data and valid transition indices from EMG_Neural_AllChannels.mat
%   2. extracts the requested emg channel
%   3. plots the emg trace in the requested time window
%   4. overlays saved valid transition points as red asterisks

% j run: plotDetectedEMGTransitions("Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData\EMG_Neural_AllChannels.mat", 1, 0, 14)

arguments
    emgNeuralFile (1,1) string
    emgChannel (1,1) double {mustBeInteger, mustBePositive}
    t0 (1,1) double
    t1 (1,1) double
end

%% ---- checks ----
if ~isfile(emgNeuralFile)
    error('file not found: %s', emgNeuralFile);
end

if t1 <= t0
    error('t1 must be greater than t0.');
end

%% ---- load file ----
S = load(emgNeuralFile);

if ~isfield(S, 'downsampEMG')
    error('downsampEMG not found in %s', emgNeuralFile);
end

if ~isfield(S, 'validTransitionsCell')
    error('validTransitionsCell not found in %s', emgNeuralFile);
end

downsampEMG = S.downsampEMG;
validTransitionsCell = S.validTransitionsCell;

%% ---- sampling rate ----
% downsampEMG in this pipeline is typically 1 khz
if isfield(S, 'fsEmg')
    fsEmg = S.fsEmg;
else
    fsEmg = 1000;
end

%% ---- get emg channel ----
% handle either channels x time or time x channels
if size(downsampEMG,1) <= 16 && size(downsampEMG,2) > size(downsampEMG,1)
    % channels x time
    if emgChannel > size(downsampEMG,1)
        error('requested emgChannel exceeds number of channels in downsampEMG.');
    end
    emgTrace = double(downsampEMG(emgChannel, :));
else
    % time x channels
    if emgChannel > size(downsampEMG,2)
        error('requested emgChannel exceeds number of channels in downsampEMG.');
    end
    emgTrace = double(downsampEMG(:, emgChannel))';
end

%% ---- get saved valid transitions for this channel ----
if emgChannel > numel(validTransitionsCell)
    error('requested emgChannel exceeds number of channels in validTransitionsCell.');
end

transitionInds = validTransitionsCell{emgChannel};
transitionInds = transitionInds(:)';
transitionInds = transitionInds(~isnan(transitionInds));
transitionInds = transitionInds(transitionInds >= 1 & transitionInds <= numel(emgTrace));

%% ---- requested time window ----
nSamp = numel(emgTrace);
timeAxis = (0:nSamp-1) ./ fsEmg;

i0 = max(1, floor(t0 * fsEmg) + 1);
i1 = min(nSamp, floor(t1 * fsEmg) + 1);

if i0 >= i1
    error('requested time window is invalid.');
end

winInds = i0:i1;
tPlot = timeAxis(winInds);
emgPlot = emgTrace(winInds);

%% ---- keep only transitions inside this window ----
inWin = transitionInds >= i0 & transitionInds <= i1;
transitionIndsWin = transitionInds(inWin);

transitionTimes = timeAxis(transitionIndsWin);
transitionVals = emgTrace(transitionIndsWin);

%% ---- plot ----
figure('Color', 'w');
plot(tPlot, emgPlot, 'LineWidth', 1.2);
hold on

if ~isempty(transitionTimes)
    plot(transitionTimes, transitionVals, 'r*', ...
        'MarkerSize', 9, ...
        'LineWidth', 1.2);
end

xlabel('Time (s)', 'FontSize', 12);
ylabel('EMG (AU)', 'FontSize', 12);

xlim([t0 t1]);

ax = gca;
ax.FontSize = 11;
ax.LineWidth = 1;
box off

%% ---- optional title ----
title(sprintf('EMG channel %d', emgChannel), 'FontSize', 13);

%% ---- command window info ----
fprintf('\nchannel %d\n', emgChannel);
fprintf('number of saved valid transitions in full trace: %d\n', numel(transitionInds));
fprintf('number of saved valid transitions in plotted window: %d\n', numel(transitionIndsWin));

if ~isempty(transitionTimes)
    fprintf('transition times in plotted window (s):\n');
    disp(transitionTimes);
end

end
