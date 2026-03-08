function plotCortexEMGandNeuralAveragesWithTransitions(dataFile, channelsToUse)
% plots two cortex-only figures:
%  1) cortex pyramidal/interneuron averages with shifted percentile bounds
%  2) full-session emg traces with detected transitions marked, plus full-session cortex firing rates

% for fig 6 of results section

% inputs
%   datafile: full path to emg_neural_allchannels.mat
%   channelstouse: vector with channel indices to pool / display (default = 1:4)

    if nargin < 2
        channelsToUse = 1:4;
    end

    % add shadederrorbar if available
    sep = 'c:\github\interneuron-analysis';
    fcn = fullfile(sep,'shadedErrorBar.m');
    if exist(fcn,'file')
        addpath(genpath(sep));
    end

    % infer base directory from the saved data file
    baseDir = fileparts(dataFile);

    % load aligned emg / neural window data and transition indices
    S = load(dataFile, ...
        'emgWindowsCell', ...
        'pyrCxWinCell', ...
        'intCxWinCell', ...
        'tAxis', ...
        'pyrCxWinShiftedCell', ...
        'intCxWinShiftedCell', ...
        'pyrCxWinShiftedMeanCell', ...
        'intCxWinShiftedMeanCell', ...
        'validTransitionsCell');

    % ---------------- 1. pool requested channels for the averaged figure ----------------
    emgPool = [];      % events x time
    pyrCxPool = [];    % events x neurons x time
    intCxPool = [];

    for ch = channelsToUse
        % keep only the triggering muscle for the emg windows from that channel
        emgChunk = squeeze(S.emgWindowsCell{ch}(:, ch, :));   % events x time
        emgPool = cat(1, emgPool, emgChunk);

        % concatenate cortex neural windows across events
        pyrCxPool = cat(1, pyrCxPool, S.pyrCxWinCell{ch});
        intCxPool = cat(1, intCxPool, S.intCxWinCell{ch});
    end

    % ---------------- 2. gather shifted cortex controls ----------------
    pyrCx_shiftMat = [];
    intCx_shiftMat = [];

    for ch = channelsToUse
        % first shift was fully saved, so average across neurons and then across events
        firstPyrCx = squeeze(mean(mean(S.pyrCxWinShiftedCell{ch}, 2, 'omitnan'), 1, 'omitnan'))';
        firstIntCx = squeeze(mean(mean(S.intCxWinShiftedCell{ch}, 2, 'omitnan'), 1, 'omitnan'))';

        % remaining 99 shifts were saved as event x time means
        pyrCxShifts = cellfun(@(M) mean(M, 1, 'omitnan'), S.pyrCxWinShiftedMeanCell(ch, :), 'UniformOutput', false);
        intCxShifts = cellfun(@(M) mean(M, 1, 'omitnan'), S.intCxWinShiftedMeanCell(ch, :), 'UniformOutput', false);

        pyrCxShifts = vertcat(pyrCxShifts{:});
        intCxShifts = vertcat(intCxShifts{:});

        pyrCxAll = [firstPyrCx; pyrCxShifts];
        intCxAll = [firstIntCx; intCxShifts];

        pyrCx_shiftMat = cat(1, pyrCx_shiftMat, pyrCxAll);
        intCx_shiftMat = cat(1, intCx_shiftMat, intCxAll);
    end

    % compute 95% shifted bounds
    pyrCx_pct = prctile(pyrCx_shiftMat, [2.5 97.5], 1);
    intCx_pct = prctile(intCx_shiftMat, [2.5 97.5], 1);

    % ---------------- 3. compute mean & sem for emg ----------------
    mEMG = mean(emgPool, 1, 'omitnan');
    seEMG = std(emgPool, 0, 1, 'omitnan') ./ sqrt(size(emgPool,1)); %#ok<NASGU>

    % ---------------- 4. compute cortex neural mean ± sem ----------------
    meanEvt = @(M) squeeze(mean(M, 2, 'omitnan'));   % average across neurons -> events x time
    grandMS = @(E) deal(mean(E, 1, 'omitnan'), std(E, 0, 1, 'omitnan') ./ sqrt(size(E,1)));

    [mPyrCx, sePyrCx] = grandMS(meanEvt(pyrCxPool));
    [mIntCx, seIntCx] = grandMS(meanEvt(intCxPool));

    % ---------------- 5. figure 1: cortex-only neural averages exactly in your yyaxis style ----------------
    figure('Name','Unshifted Cortex Neural Activity with Shift Percentile Bounds','Color','w');
    hold on;

    yyaxis left
    shadedErrorBar(S.tAxis, mPyrCx, sePyrCx, 'lineProps', {'b', 'LineWidth', 1.5});
    plot(S.tAxis, pyrCx_pct(1,:), 'b--');
    plot(S.tAxis, pyrCx_pct(2,:), 'b--');
    ylabel('Pyramidal');

    yyaxis right
    shadedErrorBar(S.tAxis, mIntCx, seIntCx, 'lineProps', {'r', 'LineWidth', 1.5});
    plot(S.tAxis, intCx_pct(1,:), 'r--');
    plot(S.tAxis, intCx_pct(2,:), 'r--');
    ylabel('Interneuron');

    title('Cortex Neural');
    xlabel('Time (ms)');
    sgtitle('unshifted Neural + 95% shifted percentile bounds');

    % ---------------- 6. load full-session emg and cortex firing rates ----------------
    E = load(fullfile(baseDir, 'EMG1ms.mat'));
    if isfield(E, 'downsampEMG')
        emgAll = E.downsampEMG;
    elseif isfield(E, 'EMG')
        emgAll = E.EMG;
    else
        error('could not find emg variable in EMG1ms.mat');
    end

    % your files are channels x time, so convert selected channels to row access
    if isvector(emgAll)
        emgAll = emgAll(:)';
    end

    F = load(fullfile(baseDir, 'NeuralFiringRates1msBins10msGauss.mat'), 'cortexFRs', 'cortexInds');
    cortexFRs = F.cortexFRs;
    cortexInds = F.cortexInds(:);

    % load centralized classifications and match this baseDir to the correct row
    conslidatedDataFoler = 'X:\David\AnalysesData';
    C = load(fullfile(conslidatedDataFoler, 'AA_classifications.mat'), 'classifications');
    classifications = C.classifications;

    animalFolders = {
        'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData', ...
        'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData', ...
        'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData', ...
        'X:\David\ArenaRecordings\D043-020425-ArenaRecording\ProcessedData'
    };

    matchRow = find(contains(string(animalFolders), string(baseDir)), 1);
    if isempty(matchRow)
        error('could not match baseDir to animalFolders for AA_classifications');
    end

    cortexLabelsAll = classifications{matchRow, 1};
    regionClass = cortexLabelsAll(cortexInds);

    intFRs = cortexFRs(regionClass == 1, :);
    pyrFRs = cortexFRs(regionClass == 0, :);

    meanIntFull = mean(intFRs, 1, 'omitnan');
    meanPyrFull = mean(pyrFRs, 1, 'omitnan');

    % collect all selected-channel transitions into one vector for xline overlays
    allTransitions = [];
    for ch = channelsToUse
        if ch <= numel(S.validTransitionsCell) && ~isempty(S.validTransitionsCell{ch})
            allTransitions = [allTransitions; S.validTransitionsCell{ch}(:)]; %#ok<AGROW>
        end
    end
    allTransitions = unique(allTransitions(:));
    allTransitions = allTransitions(~isnan(allTransitions));

    % ---------------- 7. figure 2: emg transitions + full-session cortex firing rates ----------------
    nCh = numel(channelsToUse);
    figure('Name','EMG Transitions and Cortex Firing Rates','Color','w');
    tiledlayout(nCh + 1, 1, 'TileSpacing', 'tight', 'Padding', 'compact');

    % -- one emg subplot per requested channel, with red asterisks on transitions --
    for iCh = 1:nCh
        ch = channelsToUse(iCh);
        nexttile; hold on;

        if ch > size(emgAll,1)
            title(sprintf('channel %d (missing)', ch));
            axis off;
            continue;
        end

        signal = emgAll(ch, :);
        plot(signal, 'k');

        if ch <= numel(S.validTransitionsCell) && ~isempty(S.validTransitionsCell{ch})
            validTransitions = S.validTransitionsCell{ch}(:);
            validTransitions = validTransitions(validTransitions >= 1 & validTransitions <= numel(signal));
            plot(validTransitions, signal(validTransitions), 'r*');
        else
            validTransitions = [];
        end

        title(sprintf('channel %d emg with detected transitions', ch));
        xlabel('time (ms)');
        ylabel('emg amplitude');
        box off;
    end

    % -- bottom subplot: full-session cortex pyramidal/interneuron mean firing rates --
    nexttile; hold on;
    plot(meanPyrFull, 'b', 'LineWidth', 1.5);
    plot(meanIntFull, 'r', 'LineWidth', 1.5);

    % mark all selected-channel transition times on the firing-rate plot
    for iT = 1:numel(allTransitions)
        xline(allTransitions(iT), ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
    end

    xlabel('time (ms)');
    ylabel('mean firing rate');
    title('full-session cortex firing rates with emg transition times');
    legend({'pyramidal','interneuron'}, 'Location', 'best');
    box off;

    sgtitle('Detected EMG Transitions and Cortex Population Activity');
end
