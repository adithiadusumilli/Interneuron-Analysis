function paperRandomDetectedEMGTransitionsSeparate(emgNeuralFile, emgChannel, winDur, nExamples)
% random EMG windows containing saved valid transitions - paper-formatted (like thesis fig 8A)
% vector pdf saved figs!
% run (D020):
% plotRandomDetectedEMGTransitionsSeparate("Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData\EMG_Neural_AllChannels.mat",1,14,6)

arguments
    emgNeuralFile (1,1) string
    emgChannel (1,1) double {mustBeInteger,mustBePositive}
    winDur (1,1) double {mustBePositive} = 14
    nExamples (1,1) double {mustBeInteger,mustBePositive} = 6
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
animalMatch = regexp(char(emgNeuralFile),'D\d+','match','once');
if isempty(animalMatch), animalMatch = 'Animal'; end
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

    % Unique time-range name: avoids overwriting examples from repeated runs
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
