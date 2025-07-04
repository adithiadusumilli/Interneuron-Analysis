function neuralWindowsEMGTransitions(folderPath, preSamples, postSamples)
% extracts neural activity windows around EMG transition points for all 8 EMG channels
% also separates neural activity by cell type (pyramidal vs interneurons)

% inputs:
%   folderPath   — session folder with EMG transition files and neural data
%   preSamples   — # of ms to include before each transition
%   postSamples  — # of ms to include after each transition


    for ch = 1:8
        % print which channel is being processed
        fprintf('\nProcessing EMG Channel %d...\n', ch);

        % build full path to EMG transition file for this channel
        emgEventFile = sprintf('EMG_Channel%d_Events.mat', ch);
        fullEventPath = fullfile(folderPath, emgEventFile);

        % skip this channel if no transition file found
        if ~isfile(fullEventPath)
            warning('Missing EMG event file for channel %d: %s. Skipping...', ch, emgEventFile);
            continue;
        end

        % load in EMG transition points 
        load(fullEventPath, 'validTransitions');  % transition times in ms

        % load video sync samples to align EMG and neural time frames
        load(fullfile(folderPath,'VideoSyncFrames.mat'), 'frameEMGSamples', 'frameNeuropixelSamples');

        % convert EMG transition times (20 kHz) to neural sample times (30 kHz)
        currentDir = pwd;     % save current working directory so we can return later
        cd(folderPath);       % change to the data session directory for loading local files
        neuralIndices = NeurEMGSync(validTransitions*20, frameEMGSamples, frameNeuropixelSamples, 'EMG');
        cd(currentDir);  % return to the original working directory after all processing

        % convert neural indices to 1 kHz (divide by 30 to account for downsampling)
        neuralIndices1kHz = round(neuralIndices / 30);

        % load neural firing rate data (1 ms bins)
        load(fullfile(folderPath,'NeuralFiringRates1msBins10msGauss.mat'), 'cortexFRs', 'cortexInds', 'striatumFRs', 'striatumInds');

        % load neuron classifications for pyramidal vs interneurons
        load(fullfile(folderPath,'AA_classifications.mat'), 'classifications');

        % iterate through brain regions
        for r = 1:2
            % assign region-specific variables
            if r == 1
                region = 'cortex';
                frMat = cortexFRs;
                inds = cortexInds;
                regionClass = classifications{1, 1};
            else
                region = 'striatum';
                frMat = striatumFRs;
                inds = striatumInds;
                regionClass = classifications{1, 2};
            end

            % get neuron classes for this region (0 = pyr, 1 = int)
            regionLabels = regionClass(inds);

            % separate pyramidal and interneuron indices
            pyrInds = find(regionLabels == 0);
            intInds = find(regionLabels == 1);

            % get number of events and total timepoints
            nEvents = length(neuralIndices1kHz);
            nTimepoints = size(frMat, 2);

            % preallocate storage for pyramidal neuron windows
            pyrWindows = cell(length(pyrInds), nEvents);

            % preallocate storage for interneuron windows
            intWindows = cell(length(intInds), nEvents);

            % loop through all transition events
            for e = 1:nEvents
                t = neuralIndices1kHz(e);

                % only extract if within valid range
                if t - preSamples >= 1 && t + postSamples <= nTimepoints
                    windowRange = (t - preSamples):(t + postSamples);

                    % extract pyramidal neuron windows
                    for n = 1:length(pyrInds)
                        pyrWindows{n, e} = frMat(pyrInds(n), windowRange);
                    end

                    % extract interneuron windows
                    for n = 1:length(intInds)
                        intWindows{n, e} = frMat(intInds(n), windowRange);
                    end
                end
            end

            % save results for this region and channel
            saveName = sprintf('%s_Channel%d_NeuralWindows_ByType.mat', region, ch);
            save(saveName, 'pyrWindows', 'intWindows', 'neuralIndices1kHz', 'preSamples', 'postSamples', 'pyrInds', 'intInds');

            fprintf('Saved %s pyramidal and interneuron windows for channel %d to %s\n', region, ch, saveName);
        end
    end

    % print completion message
    fprintf('\nCompleted neural window extraction by cell type for all EMG channels.\n');
end
