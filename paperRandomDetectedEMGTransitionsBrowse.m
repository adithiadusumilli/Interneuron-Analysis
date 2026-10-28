
function paperRandomDetectedEMGTransitionsBrowse(emgNeuralFile, emgChannel, winDur, nExamples, mode, selectedStart, selectedEnd)
% EMG w saved valid-transition markers: random examples + full-recording
% browse mode + manually selected pub window (run w proper input)

% orig rand mode:
% paperRandomDetectedEMGTransitionsBrowse(file,1,14,6)

% browse entire recording:
% paperRandomDetectedEMGTransitionsBrowse(file,1,14,6,"browse")
% once zoomed to a good interval, run: xlim(gca)
% record the two values (seconds) for the selected mode below.

% save exact chosen window as VECTOR PDF + FIG:
% paperRandomDetectedEMGTransitionsBrowse(file,1,14,6,"selected",120.5,128.5)

% this func plots EMG & saved detected transitions only no firing rates etc (that's the other func)

arguments
    emgNeuralFile (1,1) string
    emgChannel (1,1) double {mustBeInteger,mustBePositive}
    winDur (1,1) double {mustBePositive} = 14
    nExamples (1,1) double {mustBeInteger,mustBePositive} = 6
    mode (1,1) string = "random"
    selectedStart (1,1) double = NaN
    selectedEnd (1,1) double = NaN
end

saveDir = "X:\David\AnalysesData\InterneuronAnalyses\AA Paper Plots\emg detected transitions and neural activity aligned to emg transition";
if ~isfolder(saveDir), mkdir(saveDir); end
fontSize = 10;
axesLineWidth = 0.5;
tickLength = [0.025 0.025];

if ~isfile(emgNeuralFile), error('Could not find file: %s',emgNeuralFile); end
parentFolder = fileparts(emgNeuralFile);
emgFile = fullfile(parentFolder,'EMG1ms.mat');
if ~isfile(emgFile), error('Could not find EMG1ms.mat in %s',parentFolder); end
Strans = load(emgNeuralFile,'validTransitionsCell');
if ~isfield(Strans,'validTransitionsCell'), error('validTransitionsCell not found.'); end
if emgChannel > numel(Strans.validTransitionsCell), error('EMG channel exceeds transition cell count.'); end
transitionInds = Strans.validTransitionsCell{emgChannel};
transitionInds = transitionInds(:);
transitionInds = transitionInds(isfinite(transitionInds));
Semg = load(emgFile,'downsampEMG','fsEmg');
if ~isfield(Semg,'downsampEMG'), error('downsampEMG not found.'); end
downsampEMG = Semg.downsampEMG;
if isfield(Semg,'fsEmg') && isscalar(Semg.fsEmg) && isfinite(Semg.fsEmg) && Semg.fsEmg > 0
    fsEmg = Semg.fsEmg;
else
    fsEmg = 1000;
end
if size(downsampEMG,1) <= 16 && size(downsampEMG,2) > size(downsampEMG,1)
    if emgChannel > size(downsampEMG,1), error('EMG channel exceeds number of rows.'); end
    emgTrace = double(downsampEMG(emgChannel,:));
else
    if emgChannel > size(downsampEMG,2), error('EMG channel exceeds number of columns.'); end
    emgTrace = double(downsampEMG(:,emgChannel))';
end
transitionInds = unique(round(transitionInds));
transitionInds = transitionInds(transitionInds>=1 & transitionInds<=numel(emgTrace));
if isempty(transitionInds), error('No saved transitions for channel %d.',emgChannel); end
nSamp = numel(emgTrace);
timeAxis = (0:nSamp-1)/fsEmg;
totalDur = timeAxis(end);
if winDur >= totalDur, error('Window duration exceeds trace duration.'); end
winSamp = round(winDur*fsEmg);
if winSamp < 2, error('Window must contain at least two samples.'); end

animalMatch = regexp(char(emgNeuralFile),'D\d+','match','once');
if isempty(animalMatch), animalMatch = 'Animal'; end

switch lower(strtrim(mode))
    case "browse"
        % For multi-hour traces, display a min/max envelope so spikes are
        % visible and MATLAB can pan/zoom without plotting millions of points.
        % Selected-window mode always plots the original full-resolution data.
        maxDisplayBins = 120000;
        blockSize = max(1,ceil(nSamp/maxDisplayBins));
        if blockSize == 1
            tBrowse = timeAxis;
            yBrowse = emgTrace;
        else
            nBlocks = ceil(nSamp/blockSize);
            tBrowse = nan(1,2*nBlocks);
            yBrowse = nan(1,2*nBlocks);
            for b = 1:nBlocks
                i0 = (b-1)*blockSize+1;
                i1 = min(b*blockSize,nSamp);
                segment = emgTrace(i0:i1);
                finiteSeg = segment(isfinite(segment));
                if isempty(finiteSeg), continue; end
                tMid = (i0+i1-2)/(2*fsEmg);
                tBrowse(2*b-1:2*b) = tMid;
                yBrowse(2*b-1:2*b) = [min(finiteSeg),max(finiteSeg)];
            end
        end
        fig = figure('Color','w', ...
            'Name',sprintf('%s EMG ch%d FULL RECORDING - BROWSE',animalMatch,emgChannel), ...
            'Position',[80 100 1400 520]);
        ax = axes(fig); hold(ax,'on');
        plot(ax,tBrowse,yBrowse,'LineWidth',0.75,'Color',[0 0.447 0.741]);
        % Red stars are plotted at a constant display height, as in the original.
        finiteTrace = emgTrace(isfinite(emgTrace));
        if isempty(finiteTrace), error('EMG trace contains no finite samples.'); end
        yStar = 0.35*max(finiteTrace);
        transitionTimes = timeAxis(transitionInds);
        plot(ax,transitionTimes,yStar*ones(size(transitionTimes)), ...
            'r*','MarkerSize',5,'LineWidth',0.8);
        xlabel(ax,'Time (seconds)','FontSize',fontSize);
        ylabel(ax,'EMG (AU)','FontSize',fontSize);
        title(ax,sprintf('%s EMG channel %d: full recording (browse with pan/zoom)', ...
            animalMatch,emgChannel),'FontSize',fontSize);
        xlim(ax,[timeAxis(1),timeAxis(end)]);
        set(ax,'FontSize',fontSize,'LineWidth',axesLineWidth,'TickDir','out', ...
            'TickLength',tickLength,'XColor','k','YColor','k','Color','w');
        grid(ax,'off'); box(ax,'off');
        % Enable MATLAB's standard axes toolbar (zoom/pan). No PDF saved here.
        fprintf('\nBROWSE MODE: Use the figure toolbar to zoom and pan.\n');
        fprintf('When you find a clean interval with 2-3 peaks, click its axes and run:\n');
        fprintf('    xlim(gca)\n');
        fprintf('Then rerun with mode "selected" and those two limits (seconds).\n');
        fprintf('Browse view is downsampled as a min/max envelope when needed;\n');
        fprintf('selected mode uses every original EMG sample.\n\n');
        return;

    case "selected"
        if ~isfinite(selectedStart) || ~isfinite(selectedEnd) || ...
                selectedStart < 0 || selectedEnd > totalDur || selectedEnd <= selectedStart
            error('Selected mode requires 0 <= selectedStart < selectedEnd <= %.3f seconds.',totalDur);
        end
        i0 = max(1,ceil(selectedStart*fsEmg)+1);
        i1 = min(nSamp,floor(selectedEnd*fsEmg)+1);
        if i1 <= i0, error('Selected interval is too short.'); end
        plotInds = i0:i1;
        tPlot = timeAxis(plotInds);
        emgPlot = emgTrace(plotInds);
        inWindow = transitionInds>=i0 & transitionInds<=i1;
        transitionTimes = timeAxis(transitionInds(inWindow));
        finitePlot = emgPlot(isfinite(emgPlot));
        if isempty(finitePlot), error('Selected EMG window contains no finite data.'); end
        yStar = 0.35*max(finitePlot);
        transitionVals = yStar*ones(size(transitionTimes));
        fig = figure('Color','w', ...
            'Name',sprintf('%s EMG ch%d SELECTED %.3f-%.3fs',animalMatch,emgChannel,selectedStart,selectedEnd), ...
            'Position',[100 100 1050 360]);
        ax = axes(fig); hold(ax,'on');
        plot(ax,tPlot,emgPlot,'LineWidth',1.2);
        if ~isempty(transitionTimes)
            plot(ax,transitionTimes,transitionVals,'r*','MarkerSize',9,'LineWidth',1.2);
        end
        xlabel(ax,'Time (seconds)','FontSize',fontSize);
        ylabel(ax,'EMG (AU)','FontSize',fontSize);
        title(ax,sprintf('%s EMG channel %d: %.1f to %.1f s', ...
            animalMatch,emgChannel,tPlot(1),tPlot(end)),'FontSize',fontSize);
        xlim(ax,[selectedStart selectedEnd]);
        set(ax,'FontSize',fontSize,'LineWidth',axesLineWidth,'TickDir','out', ...
            'TickLength',tickLength,'XColor','k','YColor','k','Color','w');
        grid(ax,'off'); box(ax,'off');
        stem = sprintf('%s_manual threshold crossing selected fig_ch%d_%.3fs_to_%.3fs', ...
            animalMatch,emgChannel,selectedStart,selectedEnd);
        stem = strrep(stem,'.','p');
        pdfFile = fullfile(saveDir,stem+".pdf");
        figFile = fullfile(saveDir,stem+".fig");
        version = 2;
        while isfile(pdfFile) || isfile(figFile)
            pdfFile = fullfile(saveDir,sprintf('%s_v%d.pdf',stem,version));
            figFile = fullfile(saveDir,sprintf('%s_v%d.fig',stem,version));
            version = version+1;
        end
        exportgraphics(fig,pdfFile,'ContentType','vector','BackgroundColor','white');
        savefig(fig,figFile);
        fprintf('Saved manual selection:\n%s\n%s\n',pdfFile,figFile);
        fprintf('Saved interval: %.3f to %.3f seconds; %d detected transitions.\n', ...
            selectedStart,selectedEnd,numel(transitionTimes));
        return;

    case "random"
        % Preserve original random-example selection and plotting below.
    otherwise
        error('Unknown mode "%s". Use "random", "browse", or "selected".',mode);
end

rng('shuffle');
startInds = [];
maxTries = 5000;
nTries = 0;
while numel(startInds) < nExamples && nTries < maxTries
    nTries = nTries + 1;
    centerInd = transitionInds(randi(numel(transitionInds)));
    jitterRange = round(0.35*winSamp);
    jitter = randi([-jitterRange,jitterRange],1);
    startInd = centerInd-round(winSamp/2)+jitter;
    startInd = max(1,min(startInd,nSamp-winSamp+1));
    endInd = startInd+winSamp-1;
    if ~any(transitionInds>=startInd & transitionInds<=endInd), continue; end
    if isempty(startInds) || all(abs(startInds-startInd)>round(0.5*winSamp))
        startInds(end+1) = startInd; %#ok<AGROW>
    end
end
if isempty(startInds), error('Could not find valid windows with saved transitions.'); end
startInds = sort(startInds);
for iWin = 1:numel(startInds)
    i0 = startInds(iWin);
    i1 = i0+winSamp-1;
    plotInds = i0:i1;
    tPlot = timeAxis(plotInds);
    emgPlot = emgTrace(plotInds);
    inWindow = transitionInds>=i0 & transitionInds<=i1;
    transitionTimes = timeAxis(transitionInds(inWindow));
    yStar = 0.35*max(emgPlot);
    transitionVals = yStar*ones(size(transitionTimes));

    fig = figure('Color','w','Name',sprintf('%s EMG ch%d %.3f-%.3fs',animalMatch,emgChannel,tPlot(1),tPlot(end)), ...
        'Position',[100 100 1050 360]);
    ax = axes(fig); hold(ax,'on');
    plot(ax,tPlot,emgPlot,'LineWidth',1.2);
    if ~isempty(transitionTimes)
        plot(ax,transitionTimes,transitionVals,'r*','MarkerSize',9,'LineWidth',1.2);
    end
    xlabel(ax,'Time (seconds)','FontSize',fontSize);
    ylabel(ax,'EMG (AU)','FontSize',fontSize);
    title(ax,sprintf('%s EMG channel %d: %.1f to %.1f s',animalMatch,emgChannel,tPlot(1),tPlot(end)), ...
        'FontSize',fontSize);
    xlim(ax,[tPlot(1) tPlot(end)]);
    set(ax,'FontSize',fontSize,'LineWidth',axesLineWidth,'TickDir','out', ...
        'TickLength',tickLength,'XColor','k','YColor','k','Color','w');
    grid(ax,'off'); box(ax,'off');

    stem = sprintf('%s_EMG_ch%d_%.3fs_to_%.3fs',animalMatch,emgChannel,tPlot(1),tPlot(end));
    stem = strrep(stem,'.','p');
    pdfFile = fullfile(saveDir,stem+".pdf");
    figFile = fullfile(saveDir,stem+".fig");
    version = 2;
    while isfile(pdfFile) || isfile(figFile)
        pdfFile = fullfile(saveDir,sprintf('%s_v%d.pdf',stem,version));
        figFile = fullfile(saveDir,sprintf('%s_v%d.fig',stem,version));
        version = version+1;
    end
    exportgraphics(fig,pdfFile,'ContentType','vector','BackgroundColor','white');
    savefig(fig,figFile);
    fprintf('Saved %s\n',pdfFile);
end
fprintf('\nOpened and saved %d figures for EMG channel %d (%.2f-s windows).\n',numel(startInds),emgChannel,winDur);
end
