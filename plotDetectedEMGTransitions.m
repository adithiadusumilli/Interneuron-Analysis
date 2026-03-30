function plotDetectedEMGTransitions(emgNeuralFile, emgChannel, t0, t1)
% plots emg trace with saved valid transition events overlaid as red asterisks

% inputs:
% emgNeuralFile: full path to EMG_Neural_AllChannels.mat
% emgChannel: emg channel to plot
% t0, t1: plot window in seconds

% loads: validTransitionsCell from EMG_Neural_AllChannels.mat downsampEMG from EMG1ms.mat in the same folder
%
% run: plotDetectedEMGTransitions("Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData\EMG_Neural_AllChannels.mat", 1, 0, 14)

arguments
    emgNeuralFile (1,1) string
    emgChannel (1,1) double {mustBeInteger, mustBePositive}
    t0 (1,1) double
    t1 (1,1) double
end

%% ---- checks ----
if ~isfile(emgNeuralFile)
    error('could not find file: %s', emgNeuralFile);
end

if t1 <= t0
    error('t1 must be greater than t0.');
end

%% ---- get folder and paired emg file ----
parentFolder = fileparts(emgNeuralFile);
emgFile = fullfile(parentFolder, 'EMG1ms.mat');

if ~isfile(emgFile)
    error('could not find EMG1ms.mat in folder: %s', parentFolder);
end

%% ---- load transitions ----
Strans = load(emgNeuralFile);

if ~isfield(Strans, 'validTransitionsCell')
    error('validTransitionsCell not found in %s', emgNeuralFile);
end

validTransitionsCell = Strans.validTransitionsCell;

if emgChannel > numel(validTransitionsCell)
    error('emgChannel exceeds number of channels in validTransitionsCell.');
end

transitionInds = validTransitionsCell{emgChannel};
transitionInds = transitionInds(:);
transitionInds = transitionInds(~isnan(transitionInds));

%% ---- load emg trace ----
Semg = load(emgFile);

if ~isfield(Semg, 'downsampEMG')
    error('downsampEMG not found in %s', emgFile);
end

downsampEMG = Semg.downsampEMG;

if isfield(Semg, 'fsEmg')
    fsEmg = Semg.fsEmg;
else
    fsEmg = 1000;
end

%% ---- extract requested channel ----
if size(downsampEMG,1) <= 16 && size(downsampEMG,2) > size(downsampEMG,1)
    % channels x time
    if emgChannel > size(downsampEMG,1)
        error('emgChannel exceeds number of rows in downsampEMG.');
    end
    emgTrace = double(downsampEMG(emgChannel, :));
else
    % time x channels
    if emgChannel > size(downsampEMG,2)
        error('emgChannel exceeds number of columns in downsampEMG.');
    end
    emgTrace = double(downsampEMG(:, emgChannel))';
end

%% ---- clean transition indices ----
transitionInds = round(transitionInds);
transitionInds = transitionInds(transitionInds >= 1 & transitionInds <= numel(emgTrace));

%% ---- build time axis ----
nSamp = numel(emgTrace);
timeAxis = (0:nSamp-1) ./ fsEmg;

i0 = max(1, floor(t0 * fsEmg) + 1);
i1 = min(nSamp, floor(t1 * fsEmg) + 1);

if i1 <= i0
    error('invalid plotting window after conversion to sample indices.');
end

plotInds = i0:i1;
tPlot = timeAxis(plotInds);
emgPlot = emgTrace(plotInds);

%% ---- keep only transitions in the plotted window ----
inWindow = transitionInds >= i0 & transitionInds <= i1;
transitionIndsWin = transitionInds(inWindow);

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

set(gca, 'FontSize', 11, 'LineWidth', 1);
box off

%% ---- optional console output ----
fprintf('\nplotted emg channel %d\n', emgChannel);
fprintf('number of transitions in plotted window: %d\n', numel(transitionIndsWin));

if ~isempty(transitionTimes)
    fprintf('transition times (s):\n');
    disp(transitionTimes');
end

end
