function finalNeuralWindowsEMGTransitions(folderPath, preSamples, postSamples)
% extracts neural activity windows around EMG transition points for all 8 EMG channels
% also separates neural activity by cell type (pyramidal vs interneurons)

% inputs:
%   folderPath   — session folder with EMG transition files and neural data
%   preSamples   — # of ms to include before each transition
%   postSamples  — # of ms to include after each transition

    % add shadedErrorBar path if not already available
    shadedPath = 'C:\Github\Interneuron-Analysis';  % ← change this to your actual path
    if exist(fullfile(shadedPath, 'shadedErrorBar.m'), 'file')
        addpath(genpath(shadedPath));
    else
      warning('shadedErrorBar.m not found at specified path. check the path.');
    end

    % load once (these do not change across channels)

    load(fullfile(folderPath,'VideoSyncFrames.mat'), 'frameEMGSamples', 'frameNeuropixelSamples');    % sync frames
    load(fullfile(folderPath,'NeuralFiringRates1msBins10msGauss.mat'), 'cortexFRs', 'cortexInds', 'striatumFRs', 'striatumInds'); % firing rates
    % load classifications for all animals
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
        error('could not match folderPath to an entry in animalFolders. check the list.');
    end

    % extract classification for this animal only
    classifications = classifications(matchRow, :);

    % pre-extract region×type matrices (constant for all channels)
    cortexPyr = cortexFRs(classifications{1,1}(cortexInds)==0 , :);
    cortexInt = cortexFRs(classifications{1,1}(cortexInds)==1 , :);
    striatPyr = striatumFRs(classifications{1,2}(striatumInds)==0 , :);
    striatInt = striatumFRs(classifications{1,2}(striatumInds)==1 , :);

    % container for all 8 channels (each cell holds a struct of results)
    allChannels = cell(1,8);

    % main loop over emg channels
    for ch = 1:8
        fprintf('\nProcessing EMG Channel %d...\n', ch);

        emgEventFile = sprintf('EMG_Channel%d_Events.mat', ch);
        fullEventPath = fullfile(folderPath, emgEventFile);
        if ~isfile(fullEventPath)
            warning('Missing EMG event file for channel %d: %s. Skipping...', ch, emgEventFile);
            continue;
        end

        % load in EMG transition points  (transition times in ms)
        load(fullEventPath, 'validTransitions', 'windows');

        % convert EMG transition times (20 kHz) to neural sample times (30 kHz)
        currentDir = pwd; % save current working directory
        cd(folderPath); % change to session folder
        neuralIndices = NeurEMGSync(validTransitions*20, frameEMGSamples, frameNeuropixelSamples, 'EMG');
        cd(currentDir); % return to original directory

        % convert neural indices to 1 kHz (divide by 30 to account for downsampling)
        neuralIndices1kHz = round(neuralIndices / 30);

        % allocate 3-D matrices for pyramidal/interneuron windows (events × neurons × time)
        nEvt   = length(neuralIndices1kHz);
        tAxis  = (-preSamples : postSamples); % 201-sample time axis
        nTPnts = numel(tAxis);

        validEvents = false(nEvt, 1);
        pyrCxWin   = nan(nEvt, size(cortexPyr,1),   nTPnts);
        pyrStrWin  = nan(nEvt, size(striatPyr,1),   nTPnts);
        intCxWin   = nan(nEvt, size(cortexInt,1), nTPnts);
        intStrWin  = nan(nEvt, size(striatInt,1), nTPnts);

        % loop through all transition events
        for e = 1:nEvt
            t = neuralIndices1kHz(e);

            % validate range
            if isnan(t) || t - preSamples < 1 || t + postSamples > size(cortexFRs,2)
                continue;  % skip invalid or edge-near events
            end
            
            rng = (t - preSamples):(t + postSamples);
            % fix: ensure rng contains only valid positive integers
            if any(rng < 1) || any(rng > size(cortexFRs,2))
                continue;
            end

            validEvents(e) = true;

            % fill 3-D matrices
            pyrCxWin(e,:,:)  = cortexPyr(:,  rng);
            pyrStrWin(e,:,:) = striatPyr(:, rng);
            intCxWin(e,:,:)  = cortexInt(:, rng);
            intStrWin(e,:,:) = striatInt(:, rng);
        end

        % save results for this channel into allChannels cell array
        neuralIndices1kHz = neuralIndices1kHz(validEvents);
        pyrCxWin = pyrCxWin(validEvents,:,:);
        pyrStrWin = pyrStrWin(validEvents,:,:);
        intCxWin = intCxWin(validEvents,:,:);
        intStrWin = intStrWin(validEvents,:,:);
        allChannels{ch} = struct('neuralIndices1kHz', neuralIndices1kHz, 'pyrCx',  pyrCxWin, 'pyrStr', pyrStrWin, 'intCx',  intCxWin, 'intStr', intStrWin, 'tAxis',  tAxis );

        % plot mean ± sem for cortex pyramidal vs interneuron vs EMG
        if exist('windows','var') % EMG windows stored in EMG file
            emgWin = windows; % emgWin: events × 201
            mEMG  = mean(emgWin,1,'omitnan');
            seEMG = std(emgWin,0,1,'omitnan') ./ sqrt(size(emgWin,1));

            % average pyramidal across neurons then events
            meanPyrCxEvt = squeeze(mean(pyrCxWin,2,'omitnan')); % events×time
            mPyr = mean(meanPyrCxEvt,1,'omitnan');
            sePyr = std(meanPyrCxEvt,0,1,'omitnan') ./ sqrt(size(meanPyrCxEvt,1));

            % average interneuron 
            meanIntCxEvt = squeeze(mean(intCxWin,2,'omitnan'));
            mInt = mean(meanIntCxEvt,1,'omitnan');
            seInt = std(meanIntCxEvt,0,1,'omitnan') ./ sqrt(size(meanIntCxEvt,1));

            figure('Name',sprintf('Channel %d Cortex Avg',ch));
            shadedErrorBar(tAxis, mEMG, seEMG,  'lineprops',{'k','linewidth',1.5}); hold on;
            shadedErrorBar(tAxis, mPyr, sePyr,  'lineprops',{'b','linewidth',1.5});
            shadedErrorBar(tAxis, mInt, seInt,  'lineprops',{'r','linewidth',1.5});
            legend('EMG','Pyramidal','Interneuron');
            xlabel('Time (ms)'); ylabel('Firing / EMG');
            title(sprintf('Cortex – Channel %d', ch));
        end

         % plot mean ± sem for striatum pyramidal vs interneuron vs EMG
        if exist('windows','var')
            emgWin = windows;  % events × time
            mEMG  = mean(emgWin,1,'omitnan');
            seEMG = std(emgWin,0,1,'omitnan') ./ sqrt(size(emgWin,1));
            
            % average pyramidal across neurons then events
            meanPyrStrEvt = squeeze(mean(pyrStrWin,2,'omitnan')); % [events × time]
            mPyrStr = mean(meanPyrStrEvt,1,'omitnan');
            sePyrStr = std(meanPyrStrEvt,0,1,'omitnan') ./ sqrt(size(meanPyrStrEvt,1));
           
            % average interneurons
            meanIntStrEvt = squeeze(mean(intStrWin,2,'omitnan')); % [events × time]
            mIntStr = mean(meanIntStrEvt,1,'omitnan');
            seIntStr = std(meanIntStrEvt,0,1,'omitnan') ./ sqrt(size(meanIntStrEvt,1));

            % plot
            figure('Name',sprintf('Channel %d Striatum Avg',ch));
            shadedErrorBar(tAxis, mEMG, seEMG,  'lineprops',{'k','linewidth',1.5}); hold on;
            shadedErrorBar(tAxis, mPyrStr, sePyrStr,  'lineprops',{'b','linewidth',1.5});
            shadedErrorBar(tAxis, mIntStr, seIntStr,  'lineprops',{'r','linewidth',1.5});
            legend('EMG','Pyramidal','Interneuron');
            xlabel('Time (ms)'); ylabel('Firing / EMG');
            title(sprintf('Striatum – Channel %d', ch));
        end
    end

    % save one combined .mat file for all channels

    save(fullfile(folderPath,'AllChannels_NeuralWindows_ByType.mat'), 'allChannels', 'preSamples', 'postSamples');

    fprintf('\nCompleted neural window extraction by cell type for all EMG channels.\n');
end
