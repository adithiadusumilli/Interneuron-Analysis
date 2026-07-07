function extractEMGandNeuralWindows_getMouseDataNames(mouseID, baseSessionName, probeRegion, threshold, baselineDur, minSeparation, preSamples, postSamples)
% combines EMG event extraction and neural window extraction around EMG transitions
% also separates neural activity by cell type (pyramidal vs interneurons)

% updated version:
%   - uses David's getMouseDataNames so D050/D054 and future animals with weird naming work
%   - uses the new AA_classifications.mat row metadata to match mouse/session to classifications
%   - keeps the EMG thresholding, shifting, neural sync, window extraction, and save variables the same

% inputs:
%   mouseID  — e.g. 'D026', 'D050'
%   baseSessionName — e.g. 'D026-032923-ArenaRecording'
%   probeRegion — '' for one-probe animals, 'CFA' or 'tjMC' for two-probe animals
%   threshold — activation threshold value (e.g., 100)
%   baselineDur — minimum quiet period (in ms) before activation (default: 500)
%   minSeparation — minimum ms between separate events (default: 500)
%   preSamples — # of ms to include before each transition (default: 500)
%   postSamples — # of ms to include after each transition (default: 500)

    if nargin < 5 || isempty(baselineDur), baselineDur = 500; end
    if nargin < 6 || isempty(minSeparation), minSeparation = 500; end
    if nargin < 7 || isempty(preSamples), preSamples = 500; end
    if nargin < 8 || isempty(postSamples), postSamples = 500; end

    %% ---------------- get file paths from David's helper ----------------
    addpath('C:\Users\mirilab\Documents\GlobusTransfer');

    dataNames = getMouseDataNames(mouseID, baseSessionName, probeRegion);
    folderPath = dataNames.processedDataFolder;

    fprintf('\nextracting EMG/neural windows for %s | %s | probeRegion=%s\n', ...
        mouseID, baseSessionName, string(probeRegion));
    fprintf('processedDataFolder: %s\n', folderPath);
    fprintf('firing rate file:    %s\n', dataNames.NeuralFiringRates1msBins10msGauss);

    %% ---------------- load EMG ----------------
    load(dataNames.EMG1ms, 'downsampEMG');  % must contain 8 x n matrix
    [nChannels, nPoints] = size(downsampEMG);
    nChannels = 4;
    fs = 1000; %#ok<NASGU> % 1 kHz sampling rate

    %% ---------------- identify good/bad EMG channels ----------------
    % badEMGChans is saved in UMAP.mat. I only care about channels 1:4 here.
    % The extraction below still runs exactly the same across all 4 channels.
    % At the very end, I save the usual variable names using only good channels,
    % and save the original all-channel versions separately as GoodAndBad variables.
    if isfile(dataNames.UMAPFile)
        U = load(dataNames.UMAPFile, 'badEMGChans');

        if isfield(U, 'badEMGChans') && ~isempty(U.badEMGChans)
            badEMGChans = U.badEMGChans(:)';
        else
            badEMGChans = [];
        end
    else
        warning('UMAP.mat not found, assuming no bad EMG channels.');
        badEMGChans = [];
    end

    badEMGChans = intersect(badEMGChans, 1:nChannels);
    goodEMGChans = setdiff(1:nChannels, badEMGChans);

    if isempty(goodEMGChans)
        error('All EMG channels 1:4 are marked bad. Cannot save good-channel-only windows.');
    end

    fprintf('bad EMG channels among 1:4: ');
    disp(badEMGChans);
    fprintf('good EMG channels used for main saved variables: ');
    disp(goodEMGChans);


    % 8-cell holders for EMG and neural results
    validTransitionsCell = cell(nChannels,1); % {ch} = [events]
    emgWindowsCell = cell(nChannels,1); % {ch} = events x 8 x timePts
    tAxis = (-preSamples:postSamples); % common time vector
    validTransitionsNeurCell = cell(nChannels,1); % {ch} = synced neural indices
    pyrCxWinCell = cell(nChannels,1); % {ch} = events x neurons x time
    intCxWinCell = cell(nChannels,1);
    pyrStrWinCell = cell(nChannels,1);
    intStrWinCell = cell(nChannels,1);

    %% ---------------- load shared neural data ----------------
    load(dataNames.VideoSyncFrames, 'frameEMGSamples', 'frameNeuropixelSamples');
    load(dataNames.NeuralFiringRates1msBins10msGauss, ...
        'cortexFRs', 'cortexInds', 'striatumFRs', 'striatumInds');

    consolidatedDataFolder = 'X:\David\AnalysesData';
    classFile = fullfile(consolidatedDataFolder, 'AA_classifications.mat');

    % Load metadata if it exists. The new AA_classifications.mat should have
    % mouseIDs/baseSessionNames saved in the same row order as classifications.
    C = load(classFile, 'classifications', 'mouseIDs', 'baseSessionNames', 'probeRegions');
    classificationsAll = C.classifications;

    if isfield(C, 'mouseIDs') && isfield(C, 'baseSessionNames')
        matchRow = find(strcmp(C.mouseIDs, mouseID) & strcmp(C.baseSessionNames, baseSessionName), 1);
    else
        % fallback for older classification files that don't have metadata saved
        animalOrder = {'D026','D020','D024','D043','D050','D054'};
        matchRow = find(strcmp(animalOrder, mouseID), 1);
    end

    if isempty(matchRow)
        error('Could not match %s / %s to a row in AA_classifications.mat.', mouseID, baseSessionName);
    end

    classifications = classificationsAll(matchRow, :); % classification row for this animal/session only

    fprintf('matched AA_classifications row: %d\n', matchRow);

    %% ---------------- pre-extract region x type matrices ----------------
    % New classifications uses NaN for both wrong-region and unclassified gap neurons,
    % so ==0 and ==1 automatically use only confidently classified cells.
    cortexPyr = cortexFRs(classifications{1,1}(cortexInds)==0 , :);
    cortexInt = cortexFRs(classifications{1,1}(cortexInds)==1 , :);
    striatPyr = striatumFRs(classifications{1,2}(striatumInds)==0 , :);
    striatInt = striatumFRs(classifications{1,2}(striatumInds)==1 , :);

    fprintf('cortex pyr=%d | cortex int=%d | striat wide=%d | striat narrow=%d\n', ...
        size(cortexPyr,1), size(cortexInt,1), size(striatPyr,1), size(striatInt,1));

    % container for all 8 channels (each cell holds a struct of results)
    allChannels = cell(1, nChannels);

    % storage for shifted EMG transition indices (8x100 cell)
    validTransitionsShiftedCell = cell(nChannels, 100); % 100 shifts

    % storing shifted neural results
    validTransitionsNeurShiftedCell = cell(nChannels, 100); % 100 shifts
    pyrCxWinShiftedCell  = cell(nChannels, 1); % only stores first shift w. full data
    intCxWinShiftedCell  = cell(nChannels, 1); % only stores first shift w. full data
    pyrStrWinShiftedCell = cell(nChannels, 1); % only stores first shift w. full data
    intStrWinShiftedCell = cell(nChannels, 1); % only stores first shift w. full data

    pyrCxWinShiftedMeanCell  = cell(nChannels, 99); % stores average of other 99 shifts
    intCxWinShiftedMeanCell  = cell(nChannels, 99); % stores average of other 99 shifts
    pyrStrWinShiftedMeanCell = cell(nChannels, 99); % stores average of other 99 shifts
    intStrWinShiftedMeanCell = cell(nChannels, 99); % stores average of other 99 shifts

    %% ---------------- main loop over EMG channels ----------------
    for ch = 1:nChannels
        tic;
        signal = downsampEMG(ch, :);  % get emg signal for current channel
        allSignal = downsampEMG(1:nChannels, :);

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
            if idx > baselineDur && all(all(allSignal(idx-baselineDur:idx-1) < threshold))
                validTransitions(end+1) = idx; %#ok<AGROW>
            end
        end
        validTransitionsCell{ch} = validTransitions;

        % extract EMG windows for both unshifted and shifted transitions
        nShifts = 100;
        minShift = 30000;  % 30 sec in ms
        emgWindowsShiftedCell = cell(nChannels, nShifts);  % 8 x 100 cell
        totalLength = nPoints;

        for s = 1:(nShifts + 1)
            if s == 1
                currentTransitions = validTransitions;
            else
                while true
                    randShift = randi([minShift, totalLength - minShift]);
                    if rand() > 0.5
                        randShift = -randShift;
                    end
                    shifted = validTransitions + randShift;

                    % wrap around: subtract totalLength + 1 if index exceeds nPoints
                    shifted(shifted > totalLength) = shifted(shifted > totalLength) - (totalLength + 1);
                    shifted(shifted < 1) = shifted(shifted < 1) + totalLength;

                    % check that they're sufficiently distant from original
                    if all(abs(shifted - validTransitions) >= minShift)
                        currentTransitions = shifted;
                        break;
                    end
                end
                validTransitionsShiftedCell{ch, s - 1} = currentTransitions;
            end

            % extract sample windows for this s
            nEvt = numel(currentTransitions);
            winLen = preSamples + postSamples + 1;
            emgWin = nan(nEvt, nChannels, winLen);   % events x 8 muscles x time

            for e = 1:nEvt
                idx = currentTransitions(e);
                if idx > preSamples && idx + postSamples <= nPoints
                    emgWin(e,:,:) = downsampEMG(1:nChannels, idx-preSamples:idx+postSamples);
                end
            end

            if s == 1
                emgWindowsCell{ch} = emgWin;
            else
                emgWindowsShiftedCell{ch, s - 1} = emgWin;
            end
        end

        % plot EMG signal and mark detected threshold crossings
        figure('Name', sprintf('%s emg channel %d', mouseID, ch));
        plot(signal); hold on;
        plot(validTransitions, signal(validTransitions), 'r*');
        title(sprintf('%s channel %d emg with detected transitions', mouseID, ch));
        xlabel('time (ms)'); ylabel('emg amplitude'); hold off;

        % ----------------------- now extracting neural windows below -----------------------
        for s = 1:(nShifts + 1)
            if s == 1
                transitions = validTransitions;
            else
                transitions = validTransitionsShiftedCell{ch, s - 1};
            end

            % sync to neural time
            currentDir = pwd;
            cd(folderPath);
            neuralIndices = NeurEMGSync(transitions * 20, frameEMGSamples, frameNeuropixelSamples, 'EMG');
            cd(currentDir);

            neuralIndices1kHz = round(neuralIndices / 30);
            validTransitionsNeurShiftedCell{ch, s} = neuralIndices1kHz;

            % allocate matrices
            nEvt   = length(neuralIndices1kHz);
            tAxis  = (-preSamples : postSamples);
            nTPnts = numel(tAxis);

            pyrCxWin   = nan(nEvt, size(cortexPyr,1),   nTPnts);
            pyrStrWin  = nan(nEvt, size(striatPyr,1),   nTPnts);
            intCxWin   = nan(nEvt, size(cortexInt,1),   nTPnts);
            intStrWin  = nan(nEvt, size(striatInt,1),   nTPnts);

            for e = 1:nEvt
                t = neuralIndices1kHz(e);
                if isnan(t) || t - preSamples < 1 || t + postSamples > size(cortexFRs,2)
                    continue;
                end
                rng = (t - preSamples):(t + postSamples);
                if any(rng < 1) || any(rng > size(cortexFRs,2)), continue; end
                pyrCxWin(e,:,:)  = cortexPyr(:,  rng);
                pyrStrWin(e,:,:) = striatPyr(:, rng);
                intCxWin(e,:,:)  = cortexInt(:, rng);
                intStrWin(e,:,:) = striatInt(:, rng);
            end

            if s == 1
                validTransitionsNeurCell{ch} = neuralIndices1kHz;
                pyrCxWinCell{ch}  = pyrCxWin;
                intCxWinCell{ch}  = intCxWin;
                pyrStrWinCell{ch} = pyrStrWin;
                intStrWinCell{ch} = intStrWin;
                allChannels{ch} = struct('neuralIndices1kHz', neuralIndices1kHz, 'pyrCx', pyrCxWin, 'pyrStr', pyrStrWin, 'intCx', intCxWin, 'intStr', intStrWin, 'tAxis', tAxis);
            else
                if s == 2  % first shift: save full matrices
                    pyrCxWinShiftedCell{ch,1}  = pyrCxWin;
                    intCxWinShiftedCell{ch,1}  = intCxWin;
                    pyrStrWinShiftedCell{ch,1} = pyrStrWin;
                    intStrWinShiftedCell{ch,1} = intStrWin;
                else
                    % for s = 3 to 101 (shifts 2-100), store only the averaged population activity
                    meanEvt = @(X) squeeze(mean(X, 2, 'omitnan'));  % avg over neurons -> events x time

                    pyrCxWinShiftedMeanCell{ch, s-2}  = meanEvt(pyrCxWin);
                    intCxWinShiftedMeanCell{ch, s-2}  = meanEvt(intCxWin);
                    pyrStrWinShiftedMeanCell{ch, s-2} = meanEvt(pyrStrWin);
                    intStrWinShiftedMeanCell{ch, s-2} = meanEvt(intStrWin);
                end
            end
        end
        disp(toc);
    end


    %% ---------------- keep all 4-channel variables before saving separate files ----------------
    % These All variables are just temporary copies in memory so I can save bad-only and good-only
    % versions with the SAME variable names in separate files. They are NOT saved as duplicate variables.
    validTransitionsCellAll = validTransitionsCell;
    validTransitionsNeurCellAll = validTransitionsNeurCell;
    validTransitionsShiftedCellAll = validTransitionsShiftedCell;
    validTransitionsNeurShiftedCellAll = validTransitionsNeurShiftedCell;

    emgWindowsCellAll = emgWindowsCell;
    emgWindowsShiftedCellAll = emgWindowsShiftedCell;

    pyrCxWinCellAll = pyrCxWinCell;
    intCxWinCellAll = intCxWinCell;
    pyrStrWinCellAll = pyrStrWinCell;
    intStrWinCellAll = intStrWinCell;

    pyrCxWinShiftedCellAll = pyrCxWinShiftedCell;
    intCxWinShiftedCellAll = intCxWinShiftedCell;
    pyrStrWinShiftedCellAll = pyrStrWinShiftedCell;
    intStrWinShiftedCellAll = intStrWinShiftedCell;

    pyrCxWinShiftedMeanCellAll = pyrCxWinShiftedMeanCell;
    intCxWinShiftedMeanCellAll = intCxWinShiftedMeanCell;
    pyrStrWinShiftedMeanCellAll = pyrStrWinShiftedMeanCell;
    intStrWinShiftedMeanCellAll = intStrWinShiftedMeanCell;

    allChannelsAll = allChannels;

    %% ---------------- save bad channels and good channels separately ----------------
    % Important: I am NOT saving duplicated GoodAndBad variables in one huge file anymore.
    % Instead:
    %   1) EMG_Neural_AllChannels_BadOnly.mat = bad channels only, same variable names
    %   2) EMG_Neural_AllChannels.mat         = good channels only, same variable names downstream uses
    %
    % This keeps the downstream code unchanged because the regular filename still contains
    % pyrCxWinCell, intCxWinCell, validTransitionsCell, etc. — just pruned to good EMG chans.
    % It also avoids the huge duplicated file that was probably causing save corruption.

    %% ---------------- save bad channels separately ----------------
    badFile = fullfile(folderPath,'EMG_Neural_AllChannels_BadOnly.mat');

    if isempty(badEMGChans)
        % no bad channels for this animal, but I still save an empty metadata file
        % so it is obvious that bad channels were checked and none existed.
        savedEMGChans = badEMGChans;

        validTransitionsCell = cell(0,1);
        validTransitionsNeurCell = cell(0,1);
        validTransitionsShiftedCell = cell(0,100);
        validTransitionsNeurShiftedCell = cell(0,100);

        emgWindowsCell = cell(0,1);
        emgWindowsShiftedCell = cell(0,100);

        pyrCxWinCell = cell(0,1);
        intCxWinCell = cell(0,1);
        pyrStrWinCell = cell(0,1);
        intStrWinCell = cell(0,1);

        pyrCxWinShiftedCell = cell(0,1);
        intCxWinShiftedCell = cell(0,1);
        pyrStrWinShiftedCell = cell(0,1);
        intStrWinShiftedCell = cell(0,1);

        pyrCxWinShiftedMeanCell = cell(0,99);
        intCxWinShiftedMeanCell = cell(0,99);
        pyrStrWinShiftedMeanCell = cell(0,99);
        intStrWinShiftedMeanCell = cell(0,99);

        allChannels = cell(1,0);

    else
        % make bad-channel-only versions using the original all-channel variables
        savedEMGChans = badEMGChans;

        validTransitionsCell = validTransitionsCellAll(badEMGChans);
        validTransitionsNeurCell = validTransitionsNeurCellAll(badEMGChans);
        validTransitionsShiftedCell = validTransitionsShiftedCellAll(badEMGChans, :);
        validTransitionsNeurShiftedCell = validTransitionsNeurShiftedCellAll(badEMGChans, :);

        emgWindowsCell = emgWindowsCellAll(badEMGChans);
        emgWindowsShiftedCell = emgWindowsShiftedCellAll(badEMGChans, :);

        pyrCxWinCell = pyrCxWinCellAll(badEMGChans);
        intCxWinCell = intCxWinCellAll(badEMGChans);
        pyrStrWinCell = pyrStrWinCellAll(badEMGChans);
        intStrWinCell = intStrWinCellAll(badEMGChans);

        pyrCxWinShiftedCell = pyrCxWinShiftedCellAll(badEMGChans, :);
        intCxWinShiftedCell = intCxWinShiftedCellAll(badEMGChans, :);
        pyrStrWinShiftedCell = pyrStrWinShiftedCellAll(badEMGChans, :);
        intStrWinShiftedCell = intStrWinShiftedCellAll(badEMGChans, :);

        pyrCxWinShiftedMeanCell = pyrCxWinShiftedMeanCellAll(badEMGChans, :);
        intCxWinShiftedMeanCell = intCxWinShiftedMeanCellAll(badEMGChans, :);
        pyrStrWinShiftedMeanCell = pyrStrWinShiftedMeanCellAll(badEMGChans, :);
        intStrWinShiftedMeanCell = intStrWinShiftedMeanCellAll(badEMGChans, :);

        allChannels = allChannelsAll(badEMGChans);
    end

    save(badFile, ...
        'validTransitionsCell', 'validTransitionsNeurCell', ...
        'validTransitionsShiftedCell', 'validTransitionsNeurShiftedCell', ...
        'emgWindowsCell', 'emgWindowsShiftedCell', ...
        'pyrCxWinCell', 'intCxWinCell', 'pyrStrWinCell', 'intStrWinCell', ...
        'pyrCxWinShiftedCell', 'intCxWinShiftedCell', 'pyrStrWinShiftedCell', 'intStrWinShiftedCell', ...
        'pyrCxWinShiftedMeanCell', 'intCxWinShiftedMeanCell', 'pyrStrWinShiftedMeanCell', 'intStrWinShiftedMeanCell', ...
        'allChannels', ...
        'badEMGChans', 'goodEMGChans', 'savedEMGChans', ...
        'preSamples', 'postSamples', 'threshold', 'baselineDur', 'minSeparation', 'tAxis', ...
        'mouseID', 'baseSessionName', 'probeRegion', 'folderPath', 'matchRow', ...
        '-v7.3');

    fprintf('\nsaved bad-channel-only file to:\n%s\n', badFile);

    %% ---------------- save good channels to regular downstream filename ----------------
    goodFile = fullfile(folderPath,'EMG_Neural_AllChannels.mat');

    % make good-channel-only versions using the original all-channel variables
    savedEMGChans = goodEMGChans;

    validTransitionsCell = validTransitionsCellAll(goodEMGChans);
    validTransitionsNeurCell = validTransitionsNeurCellAll(goodEMGChans);
    validTransitionsShiftedCell = validTransitionsShiftedCellAll(goodEMGChans, :);
    validTransitionsNeurShiftedCell = validTransitionsNeurShiftedCellAll(goodEMGChans, :);

    emgWindowsCell = emgWindowsCellAll(goodEMGChans);
    emgWindowsShiftedCell = emgWindowsShiftedCellAll(goodEMGChans, :);

    pyrCxWinCell = pyrCxWinCellAll(goodEMGChans);
    intCxWinCell = intCxWinCellAll(goodEMGChans);
    pyrStrWinCell = pyrStrWinCellAll(goodEMGChans);
    intStrWinCell = intStrWinCellAll(goodEMGChans);

    pyrCxWinShiftedCell = pyrCxWinShiftedCellAll(goodEMGChans, :);
    intCxWinShiftedCell = intCxWinShiftedCellAll(goodEMGChans, :);
    pyrStrWinShiftedCell = pyrStrWinShiftedCellAll(goodEMGChans, :);
    intStrWinShiftedCell = intStrWinShiftedCellAll(goodEMGChans, :);

    pyrCxWinShiftedMeanCell = pyrCxWinShiftedMeanCellAll(goodEMGChans, :);
    intCxWinShiftedMeanCell = intCxWinShiftedMeanCellAll(goodEMGChans, :);
    pyrStrWinShiftedMeanCell = pyrStrWinShiftedMeanCellAll(goodEMGChans, :);
    intStrWinShiftedMeanCell = intStrWinShiftedMeanCellAll(goodEMGChans, :);

    allChannels = allChannelsAll(goodEMGChans);

    save(goodFile, ...
        'validTransitionsCell', 'validTransitionsNeurCell', ...
        'validTransitionsShiftedCell', 'validTransitionsNeurShiftedCell', ...
        'emgWindowsCell', 'emgWindowsShiftedCell', ...
        'pyrCxWinCell', 'intCxWinCell', 'pyrStrWinCell', 'intStrWinCell', ...
        'pyrCxWinShiftedCell', 'intCxWinShiftedCell', 'pyrStrWinShiftedCell', 'intStrWinShiftedCell', ...
        'pyrCxWinShiftedMeanCell', 'intCxWinShiftedMeanCell', 'pyrStrWinShiftedMeanCell', 'intStrWinShiftedMeanCell', ...
        'allChannels', ...
        'badEMGChans', 'goodEMGChans', 'savedEMGChans', ...
        'preSamples', 'postSamples', 'threshold', 'baselineDur', 'minSeparation', 'tAxis', ...
        'mouseID', 'baseSessionName', 'probeRegion', 'folderPath', 'matchRow', ...
        '-v7.3');

    fprintf('\nsaved good-channel-only downstream file to:\n%s\n', goodFile);

end
