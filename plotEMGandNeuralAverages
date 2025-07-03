function plotEMGandNeuralAverages(dataFile, channelsToUse)
%  plots mean ± sem traces for emg, cortex-pyramidal, cortex-interneuron, striatum-pyramidal, and striatum-interneuron
%  data come from the output of extractEMGandNeuralWindows
%  only events originating from the specified emg channels are pooled.

%  inputs
%    datafile       : full path to 'emg_neural_allchannels.mat'
%    channelstouse  : vector with channel indices to pool (default = 1:4)

    if nargin < 2, channelsToUse = 1:4; end          % default: use channels 1-4

    % add shadederrorbar (edit the path if you keep it elsewhere)
    sep = 'c:\github\interneuron-analysis';
    fcn = fullfile(sep,'shadedErrorBar.m');
    if exist(fcn,'file'), addpath(genpath(sep)); end

    % load the variables saved by extractEMGandNeuralWindows
    S = load(dataFile, 'emgWindowsCell','pyrCxWinCell','intCxWinCell', 'pyrStrWinCell','intStrWinCell','tAxis', 'pyrCxWinShiftedCell', 'intCxWinShiftedCell', 'pyrStrWinShiftedCell', 'intStrWinShiftedCell', 'pyrCxWinShiftedMeanCell','intCxWinShiftedMeanCell','pyrStrWinShiftedMeanCell','intStrWinShiftedMeanCell');

    % -------------- 1. pool events from the requested emg channels --------------
    emgPool = [];   % events × time
    pyrCxPool = [];   % events × neurons × time
    intCxPool = [];
    pyrStrPool = [];
    intStrPool = [];

    for ch = channelsToUse
        % grab emg windows for this channel → keep only the triggering muscle (dim-2 = ch)
        emgChunk = squeeze(S.emgWindowsCell{ch}(:, ch, :));    % (events × time)
        emgPool  = cat(1, emgPool, emgChunk);                    % concatenate along events

        % concatenate neural windows along the event dimension
        pyrCxPool   = cat(1 , pyrCxPool  , S.pyrCxWinCell{ch});  % (events × neurons × time)
        intCxPool   = cat(1 , intCxPool  , S.intCxWinCell{ch});
        pyrStrPool  = cat(1 , pyrStrPool , S.pyrStrWinCell{ch});
        intStrPool  = cat(1 , intStrPool , S.intStrWinCell{ch});
    end

    % --- control condition: Shifted neural (s=2), unshifted EMG ---
    % --- collect all 100 shifted means per channel (1 full + 99 stored) ---
    pyrCx_shiftMat = [];
    intCx_shiftMat = [];
    pyrStr_shiftMat = [];
    intStr_shiftMat = [];

    for ch = channelsToUse
        % first shift: mean across events & # neurons → 1 × time
        firstPyrCx  = squeeze(mean(mean(S.pyrCxWinShiftedCell{ch}, 2, 'omitnan'), 1, 'omitnan'))';  % → 1×time
        firstIntCx  = squeeze(mean(mean(S.intCxWinShiftedCell{ch}, 2, 'omitnan'), 1, 'omitnan'))';
        firstPyrStr = squeeze(mean(mean(S.pyrStrWinShiftedCell{ch}, 2, 'omitnan'), 1, 'omitnan'))';
        firstIntStr = squeeze(mean(mean(S.intStrWinShiftedCell{ch}, 2, 'omitnan'), 1, 'omitnan'))';

        % remaining 99 shifts (stored as 1x99 cells of 1x1001 each)
        pyrCxShifts  = S.pyrCxWinShiftedMeanCell(ch, :);
        intCxShifts  = S.intCxWinShiftedMeanCell(ch, :);
        pyrStrShifts = S.pyrStrWinShiftedMeanCell(ch, :);
        intStrShifts = S.intStrWinShiftedMeanCell(ch, :);

        % concatenate into 1x100 cell array
        pyrCxAll = [ {firstPyrCx},  pyrCxShifts ];
        intCxAll = [ {firstIntCx},  intCxShifts ];
        pyrStrAll = [ {firstPyrStr}, pyrStrShifts ];
        intStrAll = [ {firstIntStr}, intStrShifts ];

        % concatenate across channels (each row = one shift)
        pyrCx_shiftMat = cat(1, pyrCx_shiftMat, vertcat(pyrCxAll{:}));
        intCx_shiftMat  = cat(1, intCx_shiftMat , vertcat(intCxAll{:}));
        pyrStr_shiftMat = cat(1, pyrStr_shiftMat, vertcat(pyrStrAll{:}));
        intStr_shiftMat = cat(1, intStr_shiftMat, vertcat(intStrAll{:}));
    end

    % --- compute 95% bounds ---
    pyrCx_pct = prctile(pyrCx_shiftMat, [2.5 97.5], 1);
    intCx_pct = prctile(intCx_shiftMat, [2.5 97.5], 1);
    pyrStr_pct = prctile(pyrStr_shiftMat, [2.5 97.5], 1);
    intStr_pct = prctile(intStr_shiftMat, [2.5 97.5], 1);

    % -------------- 2. compute mean & sem for emg  (events × time) --------------
    mEMG  = mean(emgPool , 1 , 'omitnan');
    seEMG = std (emgPool , 0 , 1 , 'omitnan') ./ sqrt(size(emgPool,1));

    % -------------- 3. compute mean & sem for neural traces --------------
    %    – first average across neurons for each event,
    %      then compute grand mean/sem across those event-means
    meanEvt = @(M) squeeze(mean(M, 2 , 'omitnan')); % helper: avg over neurons → (events × time)
    grandMS = @(E) deal(mean(E,1,'omitnan'), std (E,0,1,'omitnan')./sqrt(size(E,1)));

    [mPyrCx , sePyrCx] = grandMS(meanEvt(pyrCxPool));
    [mIntCx , seIntCx] = grandMS(meanEvt(intCxPool));
    [mPyrStr, sePyrStr] = grandMS(meanEvt(pyrStrPool));
    [mIntStr, seIntStr] = grandMS(meanEvt(intStrPool));

    % -------------- 4. plot with tiledlayout  (2 rows × 3 columns) --------------
    % unified cortex/striatum plots with shifted percentiles
    figure('Name','Unshifted Neural Activity with Shift Percentile Bounds','Color','w');
    tl = tiledlayout(2, 2, 'TileSpacing', 'tight', 'Padding', 'compact');

    % -- cortex neural --
    nexttile(1); hold on;
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

    % -- striatum neural --
    nexttile(2); hold on;
    yyaxis left
    shadedErrorBar(S.tAxis, mPyrStr, sePyrStr, 'lineProps', {'b', 'LineWidth', 1.5});
    plot(S.tAxis, pyrStr_pct(1,:), 'b--');
    plot(S.tAxis, pyrStr_pct(2,:), 'b--');
    ylabel('Pyramidal');

    yyaxis right
    shadedErrorBar(S.tAxis, mIntStr, seIntStr, 'lineProps', {'r', 'LineWidth', 1.5});
    plot(S.tAxis, intStr_pct(1,:), 'r--');
    plot(S.tAxis, intStr_pct(2,:), 'r--');
    ylabel('Interneuron');

    title('Striatum Neural');
    xlabel('Time (ms)');

    sgtitle('unshifted Neural + 95% shifted percentile bounds');
end
