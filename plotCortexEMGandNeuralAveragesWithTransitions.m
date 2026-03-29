function plotCortexEMGandNeuralAveragesWithTransitions(dataFile, channelsToUse, doBaselineNorm)
% plots three cortex-only figures:
%  1) cortex pyramidal/interneuron averages aligned to emg transition events with shifted percentile bounds
%  2) zoomed channel 1 emg transition panel + cortex firing rates
%  3) cortex neural activity aligned to emg transitions with mean emg subplot underneath

% inputs
%   datafile: full path to emg_neural_allchannels.mat
%   channelstouse: vector with channel indices to pool for the aligned averages (default = 1:4)
%   doBaselineNorm: logical, per-trial baseline subtract (-500 to -450 ms)
%                   after avg across neurons but before avg across trials
%                   (default = true)

    if nargin < 2
        channelsToUse = 1:4;
    end
    if nargin < 3 || isempty(doBaselineNorm)
        doBaselineNorm = true;
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

    % ---------------- helper funcs ----------------
    meanEvt = @(M) squeeze(mean(M, 2, 'omitnan'));
    grandMS = @(E) deal(mean(E, 1, 'omitnan'), std(E, 0, 1, 'omitnan') ./ sqrt(size(E,1)));

    function Eout = subtractTrialBaseline(Ein, tAxisLocal, tStart, tEnd)
        [~, iStart] = min(abs(tAxisLocal - tStart));
        [~, iEnd]   = min(abs(tAxisLocal - tEnd));
        if iEnd < iStart
            tmp = iStart; iStart = iEnd; iEnd = tmp;
        end
        baselines = mean(Ein(:, iStart:iEnd), 2, 'omitnan');
        Eout = Ein - baselines;
    end

    % ---------------- 1. pool requested channels for the averaged figures ----------------
    emgPool = [];
    pyrCxPool = [];
    intCxPool = [];

    for ch = channelsToUse
        emgChunk = squeeze(S.emgWindowsCell{ch}(:, ch, :));
        emgPool = cat(1, emgPool, emgChunk);

        pyrCxPool = cat(1, pyrCxPool, S.pyrCxWinCell{ch});
        intCxPool = cat(1, intCxPool, S.intCxWinCell{ch});
    end

    % ---------------- 2. gather shifted cortex controls ----------------
    pyrCx_shiftMat = [];
    intCx_shiftMat = [];

    for ch = channelsToUse
        firstPyrCx_evt = meanEvt(S.pyrCxWinShiftedCell{ch});
        firstIntCx_evt = meanEvt(S.intCxWinShiftedCell{ch});

        if doBaselineNorm
            firstPyrCx_evt = subtractTrialBaseline(firstPyrCx_evt, S.tAxis, -500, -450);
            firstIntCx_evt = subtractTrialBaseline(firstIntCx_evt, S.tAxis, -500, -450);
        end

        firstPyrCx = mean(firstPyrCx_evt, 1, 'omitnan');
        firstIntCx = mean(firstIntCx_evt, 1, 'omitnan');

        pyrCxShifts = cell(size(S.pyrCxWinShiftedMeanCell, 2), 1);
        intCxShifts = cell(size(S.intCxWinShiftedMeanCell, 2), 1);

        for k = 1:size(S.pyrCxWinShiftedMeanCell, 2)
            E = S.pyrCxWinShiftedMeanCell{ch, k};
            if doBaselineNorm
                E = subtractTrialBaseline(E, S.tAxis, -500, -450);
            end
            pyrCxShifts{k} = mean(E, 1, 'omitnan');
        end

        for k = 1:size(S.intCxWinShiftedMeanCell, 2)
            E = S.intCxWinShiftedMeanCell{ch, k};
            if doBaselineNorm
                E = subtractTrialBaseline(E, S.tAxis, -500, -450);
            end
            intCxShifts{k} = mean(E, 1, 'omitnan');
        end

        pyrCxShifts = vertcat(pyrCxShifts{:});
        intCxShifts = vertcat(intCxShifts{:});

        pyrCxAll = [firstPyrCx; pyrCxShifts];
        intCxAll = [firstIntCx; intCxShifts];

        pyrCx_shiftMat = cat(1, pyrCx_shiftMat, pyrCxAll);
        intCx_shiftMat = cat(1, intCx_shiftMat, intCxAll);
    end

    pyrCx_pct = prctile(pyrCx_shiftMat, [2.5 97.5], 1);
    intCx_pct = prctile(intCx_shiftMat, [2.5 97.5], 1);

    % ---------------- 3. compute mean & sem for emg ----------------
    mEMG = mean(emgPool, 1, 'omitnan');
    seEMG = std(emgPool, 0, 1, 'omitnan') ./ sqrt(size(emgPool,1));

    % ---------------- 4. compute cortex neural mean ± sem ----------------
    pyrCx_evt = meanEvt(pyrCxPool);
    intCx_evt = meanEvt(intCxPool);

    if doBaselineNorm
        pyrCx_evt = subtractTrialBaseline(pyrCx_evt, S.tAxis, -500, -450);
        intCx_evt = subtractTrialBaseline(intCx_evt, S.tAxis, -500, -450);
    end

    [mPyrCx, sePyrCx] = grandMS(pyrCx_evt);
    [mIntCx, seIntCx] = grandMS(intCx_evt);

    % ---------------- 5. figure 1 ----------------
    figure('Name','Cortex Neural Activity Aligned to EMG Transition Events','Color','w');
    hold on;

    yyaxis left
    hPyr = shadedErrorBar(S.tAxis, mPyrCx, sePyrCx, 'lineProps', {'b', 'LineWidth', 1.5});
    hPyrLo = plot(S.tAxis, pyrCx_pct(1,:), 'b--', 'LineWidth', 1.0);
    hPyrHi = plot(S.tAxis, pyrCx_pct(2,:), 'b--', 'LineWidth', 1.0);
    ylabel('Pyramidal Firing Rate', 'FontSize', 18);

    yyaxis right
    hInt = shadedErrorBar(S.tAxis, mIntCx, seIntCx, 'lineProps', {'r', 'LineWidth', 1.5});
    hIntLo = plot(S.tAxis, intCx_pct(1,:), 'r--', 'LineWidth', 1.0);
    hIntHi = plot(S.tAxis, intCx_pct(2,:), 'r--', 'LineWidth', 1.0);
    ylabel('Interneuron Firing Rate', 'FontSize', 18);

    xlabel('Time Relative to EMG Transition (ms)', 'FontSize', 18);

    ax1 = gca;
    ax1.FontSize = 16;
    ax1.LineWidth = 1;
    ax1.TickDir = 'out';
    box off;

    yyaxis left
    ax1.YAxis(1).Exponent = 0;
    yyaxis right
    ax1.YAxis(2).Exponent = 0;

    t1 = title('Cortex Neural Activity Aligned to EMG Transition Events', 'FontSize', 18);
    t1.Units = 'normalized';
    t1.Position(2) = 0.98;

    legend([hPyr.mainLine, hPyrLo, hInt.mainLine, hIntLo], ...
        {'Pyramidal Mean \pm SEM', ...
         'Pyramidal Shifted 95% Bounds', ...
         'Interneuron Mean \pm SEM', ...
         'Interneuron Shifted 95% Bounds'}, ...
        'Location', 'best');

    % ---------------- 6. load full-session emg and cortex firing rates ----------------
    E = load(fullfile(baseDir, 'EMG1ms.mat'));
    if isfield(E, 'downsampEMG')
        emgAll = E.downsampEMG;
    elseif isfield(E, 'EMG')
        emgAll = E.EMG;
    else
        error('could not find emg variable in EMG1ms.mat');
    end

    if isvector(emgAll)
        emgAll = emgAll(:)';
    end

    F = load(fullfile(baseDir, 'NeuralFiringRates1msBins10msGauss.mat'), 'cortexFRs', 'cortexInds');
    cortexFRs = F.cortexFRs;
    cortexInds = F.cortexInds(:);

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

    % ---------------- 7. figure 2: fixed absolute-time window with transitions ----------------
    selectedChannel = 1;

    if selectedChannel > size(emgAll,1)
        error('channel 1 does not exist in EMG data.');
    end

    signalCh1 = emgAll(selectedChannel, :);

    if selectedChannel <= numel(S.validTransitionsCell) && ~isempty(S.validTransitionsCell{selectedChannel})
        transitionsCh1 = S.validTransitionsCell{selectedChannel}(:);
        transitionsCh1 = transitionsCh1(~isnan(transitionsCh1));
        transitionsCh1 = transitionsCh1(transitionsCh1 >= 1 & transitionsCh1 <= numel(signalCh1));
        transitionsCh1 = sort(unique(transitionsCh1));
    else
        error('no valid transitions found for channel 1.');
    end

    % choose a fixed window that contains detected transitions but avoids the big spike region
    xWin = [1.5144e6 1.5150e6];
    inWin = transitionsCh1 >= xWin(1) & transitionsCh1 <= xWin(2);

    if ~any(inWin)
        error('no detected transitions found in the requested fixed window. choose a different xWin.');
    end

    figure('Name','Channel 1 EMG Transitions and Cortex Population Activity','Color','w');
    tiledlayout(2, 1, 'TileSpacing', 'tight', 'Padding', 'compact');

    % -- top subplot: channel 1 emg with detected transitions --
    nexttile; hold on;
    hEmg = plot(1:numel(signalCh1), signalCh1, 'k');

    hTrans = plot(transitionsCh1(inWin), signalCh1(transitionsCh1(inWin)), 'r*', ...
        'MarkerSize', 8, 'LineWidth', 1.2);

    title('Channel 1 EMG Detected Transitions', 'FontSize', 16);
    xlabel('Time (ms)', 'FontSize', 16);
    ylabel('Emg Amplitude', 'FontSize', 16);
    xlim(xWin);
    box off;

    ax2 = gca;
    ax2.FontSize = 14;
    ax2.LineWidth = 1;
    ax2.TickDir = 'out';

    legend([hEmg hTrans], {'EMG Signal', 'Detected Transition'}, 'Location', 'best');

    % -- bottom subplot: M1 Neural Firing Rates --
    nexttile; hold on;
    plot(1:numel(meanPyrFull), meanPyrFull, 'b', 'LineWidth', 1.5);
    plot(1:numel(meanIntFull), meanIntFull, 'r', 'LineWidth', 1.5);

    theseTransitions = transitionsCh1(inWin);
    for iT = 1:numel(theseTransitions)
        xline(theseTransitions(iT), ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.75);
    end

    xlabel('Time (ms)', 'FontSize', 16);
    ylabel('Mean Firing Rate', 'FontSize', 16);
    title('M1 Firing Rates Corresponding to EMG Channel 1 Detected Transitions', 'FontSize', 16);
    legend({'Pyramidal Neuron','Interneuron'}, 'Location', 'best');
    xlim(xWin);
    box off;

    ax3 = gca;
    ax3.FontSize = 14;
    ax3.LineWidth = 1;
    ax3.TickDir = 'out';

    % ---------------- 8. figure 3 ----------------
    figure('Name','M1 Neural Activity with EMG Subplot','Color','w');
    tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    axNeural = nexttile; hold on;

    yyaxis left
    hPyr2 = shadedErrorBar(S.tAxis, mPyrCx, sePyrCx, 'lineProps', {'b', 'LineWidth', 1.5});
    hPyrLo2 = plot(S.tAxis, pyrCx_pct(1,:), 'b--', 'LineWidth', 1.0);
    hPyrHi2 = plot(S.tAxis, pyrCx_pct(2,:), 'b--', 'LineWidth', 1.0);
    ylabel('Pyramidal Neuron Firing Rate (spikes/s)','FontSize',18)

    yyaxis right
    hInt2 = shadedErrorBar(S.tAxis, mIntCx, seIntCx, 'lineProps', {'r', 'LineWidth', 1.5});
    hIntLo2 = plot(S.tAxis, intCx_pct(1,:), 'r--', 'LineWidth', 1.0);
    hIntHi2 = plot(S.tAxis, intCx_pct(2,:), 'r--', 'LineWidth', 1.0);
    ylabel('Interneuron Firing Rate (spikes/s)','FontSize',18)

    xline(0,'k:','LineWidth',1)
    title('M1 Neural Activity During EMG Transition Windows','FontSize',18)

    axNeural.FontSize = 16;
    axNeural.LineWidth = 1;
    axNeural.TickDir = 'out';
    box off

    yyaxis left
    axNeural.YAxis(1).Exponent = 0;
    yyaxis right
    axNeural.YAxis(2).Exponent = 0;

    lgd = legend(axNeural, ...
        [hPyr2.mainLine hPyrLo2 hInt2.mainLine hIntLo2], ...
        {'Pyramidal Mean ± SEM', ...
         'Pyramidal Shifted 95% Bounds', ...
         'Interneuron Mean ± SEM', ...
         'Interneuron Shifted 95% Bounds'}, ...
        'Location','southoutside', ...
        'Orientation','horizontal');

    lgd.FontSize = 14;
    lgd.Box = 'off';

    axEMG = nexttile; hold on;

    shadedErrorBar(S.tAxis, mEMG, seEMG, 'lineProps', {'k','LineWidth',1.5});
    xline(0,'k:','LineWidth',1)

    xlabel('Time Relative to EMG Transition (ms)','FontSize',18)
    ylabel('EMG Amplitude','FontSize',18)
    title('Mean EMG Aligned to Transition Events','FontSize',18)

    axEMG.FontSize = 16;
    axEMG.LineWidth = 1;
    axEMG.TickDir = 'out';
    box off
end
