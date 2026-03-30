function plotRandomDetectedEMGTransitionsSeparate(emgNeuralFile, emgChannel, winDur, nExamples)
% plots multiple separate emg figures from random time windows that contain
% at least one saved valid transition

% inputs:
%   emgNeuralFile: full path to EMG_Neural_AllChannels.mat
%   emgChannel: emg channel to plot
%   winDur: window duration in seconds
%   nExamples: number of separate example figures to open

% func loads validTransitionsCell from EMG_Neural_AllChannels.mat & downsampEMG from EMG1ms.mat 

% each fig has 1 randomly chosen time block & overlays saved valid transitions as red asterisks

% j run: plotRandomDetectedEMGTransitionsSeparate("Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData\EMG_Neural_AllChannels.mat", 1, 14, 10)

arguments
    emgNeuralFile (1,1) string
    emgChannel (1,1) double {mustBeInteger, mustBePositive}
    winDur (1,1) double {mustBePositive} = 14
    nExamples (1,1) double {mustBeInteger, mustBePositive} = 6
end

%% ---- checks ----
if ~isfile(emgNeuralFile)
    error('could not find file: %s', emgNeuralFile);
end

%% ---- paired emg file ----
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

%% ---- load emg ----
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

%% ---- extract emg channel ----
if size(downsampEMG,1) <= 16 && size(downsampEMG,2) > size(downsampEMG,1)
    if emgChannel > size(downsampEMG,1)
        error('emgChannel exceeds number of rows in downsampEMG.');
    end
    emgTrace = double(downsampEMG(emgChannel, :));
else
    if emgChannel > size(downsampEMG,2)
        error('emgChannel exceeds number of columns in downsampEMG.');
    end
    emgTrace = double(downsampEMG(:, emgChannel))';
end

%% ---- clean transition indices ----
transitionInds = round(transitionInds);
transitionInds = transitionInds(transitionInds >= 1 & transitionInds <= numel(emgTrace));

if isempty(transitionInds)
    error('no saved transitions found for emg channel %d.', emgChannel);
end

%% ---- overall timing ----
nSamp = numel(emgTrace);
timeAxis = (0:nSamp-1) ./ fsEmg;
totalDur = timeAxis(end);

if winDur >= totalDur
    error('window duration %.2f s is longer than the total trace duration %.2f s.', winDur, totalDur);
end

winSamp = round(winDur * fsEmg);

%% ---- randomly choose windows that contain at least one transition ----
rng('shuffle');

startInds = [];
maxTries = 5000;
nFound = 0;
nTries = 0;

while nFound < nExamples && nTries < maxTries
    nTries = nTries + 1;

    centerInd = transitionInds(randi(numel(transitionInds)));

    jitterRange = round(0.35 * winSamp);
    jitter = randi([-jitterRange, jitterRange], 1);

    startInd = centerInd - round(winSamp/2) + jitter;
    startInd = max(1, startInd);
    startInd = min(startInd, nSamp - winSamp + 1);

    endInd = startInd + winSamp - 1;

    hasTransition = any(transitionInds >= startInd & transitionInds <= endInd);

    if ~hasTransition
        continue
    end

    if isempty(startInds) || all(abs(startInds - startInd) > round(0.5 * winSamp))
        startInds(end+1) = startInd; %#ok<AGROW>
        nFound = nFound + 1;
    end
end

if isempty(startInds)
    error('could not find any valid windows with saved transitions.');
end

startInds = sort(startInds);
nPlots = numel(startInds);

%% ---- make separate figures ----
for iWin = 1:nPlots
    i0 = startInds(iWin);
    i1 = i0 + winSamp - 1;

    plotInds = i0:i1;
    tPlot = timeAxis(plotInds);
    emgPlot = emgTrace(plotInds);

    inWindow = transitionInds >= i0 & transitionInds <= i1;
    transitionIndsWin = transitionInds(inWindow);
    transitionTimes = timeAxis(transitionIndsWin);

    % fixed-height asterisks like the example image
    yStar = 0.35 * max(emgPlot);
    transitionVals = yStar * ones(size(transitionTimes));

    figure('Color', 'w');
    plot(tPlot, emgPlot, 'LineWidth', 1.2);
    hold on

    if ~isempty(transitionTimes)
        plot(transitionTimes, transitionVals, 'r*', ...
            'MarkerSize', 9, ...
            'LineWidth', 1.2);
    end

    xlabel('Time (s)');
    ylabel('EMG (AU)');
    xlim([tPlot(1) tPlot(end)]);
    set(gca, 'FontSize', 11, 'LineWidth', 1);
    box off

    title(sprintf('EMG channel %d: %.1f to %.1f s', emgChannel, tPlot(1), tPlot(end)));
end

%% ---- command window output ----
fprintf('\nopened %d separate figure(s) for emg channel %d\n', nPlots, emgChannel);
fprintf('window duration: %.2f s\n', winDur);

end
