function plotDetectedEMGTransitions(baseDir, emgChannel, t0, t1)
% plots emg trace + overlays saved valid transitions (red asterisks)

% inputs:
%   baseDir: path to ProcessedData folder
%   emgChannel: emg channel to plot
%   t0, t1: time window in seconds

% uses: EMG1ms.mat for emg signal & EMG_Neural_AllChannels.mat bc it has the validTransitionsCell from 1st extract func

arguments
    baseDir (1,1) string
    emgChannel (1,1) double {mustBeInteger, mustBePositive}
    t0 (1,1) double
    t1 (1,1) double
end

%% ---- load EMG ----
emgFile = fullfile(baseDir, 'EMG1ms.mat');
S1 = load(emgFile);

if ~isfield(S1, 'downsampEMG')
    error('downsampEMG not found in EMG1ms.mat');
end

downsampEMG = S1.downsampEMG;

if isfield(S1, 'fsEmg')
    fsEmg = S1.fsEmg;
else
    fsEmg = 1000;
end

%% ---- load transitions ----
transFile = fullfile(baseDir, 'EMG_Neural_AllChannels.mat');
S2 = load(transFile);

if ~isfield(S2, 'validTransitionsCell')
    error('validTransitionsCell not found');
end

validTransitionsCell = S2.validTransitionsCell;

%% ---- extract EMG channel ----
if size(downsampEMG,1) <= 16
    emgTrace = double(downsampEMG(emgChannel, :));
else
    emgTrace = double(downsampEMG(:, emgChannel))';
end

%% ---- get transitions ----
transitionInds = validTransitionsCell{emgChannel};
transitionInds = transitionInds(~isnan(transitionInds));

%% ---- time axis ----
nSamp = numel(emgTrace);
timeAxis = (0:nSamp-1) ./ fsEmg;

%% ---- window ----
i0 = max(1, floor(t0 * fsEmg) + 1);
i1 = min(nSamp, floor(t1 * fsEmg) + 1);

tPlot = timeAxis(i0:i1);
emgPlot = emgTrace(i0:i1);

%% ---- filter transitions in window ----
inWin = transitionInds >= i0 & transitionInds <= i1;
transitionIndsWin = transitionInds(inWin);

transitionTimes = timeAxis(transitionIndsWin);
transitionVals = emgTrace(transitionIndsWin);

%% ---- plot ----
figure('Color','w');
plot(tPlot, emgPlot, 'LineWidth', 1.2);
hold on

plot(transitionTimes, transitionVals, 'r*', ...
    'MarkerSize', 9, 'LineWidth', 1.2);

xlabel('Time (s)');
ylabel('EMG (AU)');
xlim([t0 t1]);

set(gca, 'FontSize', 11, 'LineWidth', 1);
box off

title(sprintf('EMG channel %d', emgChannel));

%% ---- debug info ----
fprintf('transitions in window: %d\n', numel(transitionTimes));

end
