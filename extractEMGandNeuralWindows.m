function extractEMGandNeuralWindows(folderPath, threshold, baselineDur, minSeparation, preSamples, postSamples)
% combines EMG event extraction and neural window extraction around EMG transitions
% also separates neural activity by cell type (pyramidal vs interneurons)

% inputs:
%   folderPath     — session folder with EMG, neural, and sync data
%   threshold      — activation threshold value (e.g., 100)
%   baselineDur    — minimum quiet period (in ms) before activation (default: 200)
%   minSeparation  — minimum ms between separate events (default: 500)
%   preSamples     — # of ms to include before each transition
%   postSamples    — # of ms to include after each transition

    if nargin < 3, baselineDur = 200; end
    if nargin < 4, minSeparation = 500; end
    if nargin < 5, preSamples = 100; end
    if nargin < 6, postSamples = 100; end

    load(fullfile(folderPath,'EMG1ms'), 'downsampEMG');  % must contain 8 × n matrix
    [nChannels, nPoints] = size(downsampEMG);
    fs = 1000;  % 1 kHz sampling rate

    % 8-cell holders for EMG and neural results
    validTransitionsCell = cell(nChannels,1); % {ch} = [events]
    emgWindowsCell = cell(nChannels,1); % {ch} = events × 8 × timePts
    tAxis = (-preSamples:postSamples); % common time vector
    validTransitionsNeurCell = cell(nChannels,1); % {ch} = synced neural indices
    pyrCxWinCell = cell(nChannels,1); % {ch} = events × neurons × time
    intCxWinCell = cell(nChannels,1);
    pyrStrWinCell = cell(nChannels,1);
    intStrWinCell = cell(nChannels,1);

    % load shared neural data
    load(fullfile(folderPath,'VideoSyncFrames.mat'), 'frameEMGSamples', 'frameNeuropixelSamples');    % sync frames
    load(fullfile(folderPath,'NeuralFiringRates1msBins10msGauss.mat'), 'cortexFRs', 'cortexInds', 'striatumFRs', 'striatumInds'); % firing rates
    conslidatedDataFoler = 'X:\David\AnalysesData';
    load(fullfile(conslidatedDataFoler, 'AA_classifications.mat'), 'classifications');

    % define known base folders (match order used when AA_classifications.mat was created)
    animalFolders = {
        'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData', ...
        'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData', ...
        'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData'
    };

    % determine which animal this folderPath corresponds to
    matchRow = find(contains(animalFolders, folderPath), 1);
    if isempty(matchRow)
        error('Could not match folderPath to an entry in animalFolders. Please check the list.');
    end
    classifications = classifications(matchRow, :); % extract classification for this animal only

    % pre-extract region×type matrices (constant for all channels)
    cortexPyr = cortexFRs(classifications{1,1}(cortexInds)==0 , :);
    cortexInt = cortexFRs(classifications{1,1}(cortexInds)==1 , :);
    striatPyr = striatumFRs(classifications{1,2}(striatumInds)==0 , :);
    striatInt = striatumFRs(classifications{1,2}(striatumInds)==1 , :);

    % container for all 8 channels (each cell holds a struct of results)
    allChannels = cell(1, nChannels);

    % main loop over emg channels
    for ch = 1:nChannels
        signal = downsampEMG(ch, :);  % get emg signal for current channel

        % apply threshold to find where signal is above threshold
        above = signal > threshold;
        transitionPoints = find(diff([0, above]) == 1);

        % remove transitions that are too close together (< minSeparation)
        if numel(transitionPoints) > 1
            diffs = diff(transitionPoints);
            keepers = [true, diffs > minSeparation];
            transitionPoints = transitionPoints(keepers);
        end

        % filter out threshold crossings that do not have quiet baseline before
        validTransitions = [];
        for i = 1:numel(transitionPoints)
            idx = transitionPoints(i);
            if idx > baselineDur && all(signal(idx-baselineDur:idx-1) < threshold)
                validTransitions(end+1) = idx;
            end
        end
        validTransitionsCell{ch} = validTransitions;

        % extract sample window (before, after, plus crossing)
        nEvt   = numel(validTransitions);
        winLen = preSamples + postSamples + 1;
        emgWin = nan(nEvt, nChannels, winLen);   % events × 8 muscles × time

        for e = 1:nEvt
            idx = validTransitions(e);
            if idx > preSamples && idx + postSamples <= nPoints
                emgWin(e,:,:) = downsampEMG(:, idx-preSamples:idx+postSamples);
            end
        end
        emgWindowsCell{ch} = emgWin; % store per-channel EMG data

        % plot emg signal and mark detected threshold crossings
        figure('Name', sprintf('emg channel %d', ch));
        plot(signal); hold on;
        plot(validTransitions, signal(validTransitions), 'r*');
        title(sprintf('channel %d emg with detected transitions', ch));
        xlabel('time (ms)'); ylabel('emg amplitude'); hold off;

        % ----------------------- now extracting neural windows below -----------------------

        % convert EMG transition times (20 kHz) to neural sample times (30 kHz)
        currentDir = pwd;
        cd(folderPath);
        neuralIndices = NeurEMGSync(validTransitions*20, frameEMGSamples, frameNeuropixelSamples, 'EMG');
        cd(currentDir);

         % convert neural indices to 1 kHz (divide by 30 to account for downsampling)
        neuralIndices1kHz = round(neuralIndices / 30);
        validTransitionsNeurCell{ch} = neuralIndices1kHz;

        % allocate 3-D matrices for pyramidal/interneuron windows (events × neurons × time)
        nEvt   = length(neuralIndices1kHz);
        tAxis  = (-preSamples : postSamples);
        nTPnts = numel(tAxis);

        pyrCxWin   = nan(nEvt, size(cortexPyr,1),   nTPnts);
        pyrStrWin  = nan(nEvt, size(striatPyr,1),   nTPnts);
        intCxWin   = nan(nEvt, size(cortexInt,1),   nTPnts);
        intStrWin  = nan(nEvt, size(striatInt,1),   nTPnts);

        % loop through all transition events
        for e = 1:nEvt
            t = neuralIndices1kHz(e);

            % validate range
            if isnan(t) || t - preSamples < 1 || t + postSamples > size(cortexFRs,2)
                continue; % skip invalid or edge-near events
            end
            rng = (t - preSamples):(t + postSamples);

            % fix: ensure rng contains only valid positive integers
            if any(rng < 1) || any(rng > size(cortexFRs,2))
                continue;
            end

            % fill 3-D matrices
            pyrCxWin(e,:,:)  = cortexPyr(:,  rng);
            pyrStrWin(e,:,:) = striatPyr(:, rng);
            intCxWin(e,:,:)  = cortexInt(:, rng);
            intStrWin(e,:,:) = striatInt(:, rng);
        end

        pyrCxWinCell{ch}  = pyrCxWin;
        intCxWinCell{ch}  = intCxWin;
        pyrStrWinCell{ch} = pyrStrWin;
        intStrWinCell{ch} = intStrWin;

        % save results for this channel
        allChannels{ch} = struct('neuralIndices1kHz', neuralIndices1kHz,'pyrCx', pyrCxWin, 'pyrStr', pyrStrWin, 'intCx', intCxWin, 'intStr', intStrWin, 'tAxis', tAxis);
    end

    % save everything together
    save(fullfile(folderPath,'EMG_Neural_AllChannels.mat'), 'validTransitionsCell','validTransitionsNeurCell', 'emgWindowsCell', 'pyrCxWinCell','intCxWinCell','pyrStrWinCell','intStrWinCell', 'preSamples','postSamples','threshold','baselineDur','minSeparation','tAxis','-v7.3');
end
