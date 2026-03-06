function plotCombinedSpikeWidthGMM(baseDirs, regionToPlot)
% plots the combined spike-width distribution across multiple animals,
% fits a 2-component gmm, overlays the two gaussian components,
% and marks the intersection point between them

% inputs:
%   baseDirs : cell array of ProcessedData folder paths
%   regionToPlot : "cortex" or "striatum"  [default = "cortex"]

% j run:
% baseDirs = {
%     'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData', ...
%     'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData', ...
%     'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData', ...
%     'X:\David\ArenaRecordings\D043-020425-ArenaRecording\ProcessedData'
% };
% plotCombinedSpikeWidthGMM(baseDirs, "cortex")

arguments
    baseDirs cell
    regionToPlot (1,1) string = "cortex"
end

numFiles = numel(baseDirs);
allWidths = [];

%% ---- collect spike widths across all animals ----
for fileIndex = 1:numFiles
    baseDir = string(baseDirs{fileIndex});

    S = load(fullfile(baseDir,'neuronDataStruct.mat'),'neuronDataStruct');
    F = load(fullfile(baseDir,'NeuralFiringRates1msBins10msGauss.mat'),'cortexInds','striatumInds');

    neuronDataStruct = S.neuronDataStruct;

    switch lower(regionToPlot)
        case "cortex"
            regionInds = F.cortexInds;
        case "striatum"
            regionInds = F.striatumInds;
        otherwise
            error('regionToPlot must be "cortex" or "striatum".');
    end

    theseWidths = [];

    for i = 1:numel(regionInds)
        waveform    = neuronDataStruct(regionInds(i)).waveforms;
        biggestChan = neuronDataStruct(regionInds(i)).biggestChan;
        ap          = waveform(:,biggestChan);

        [~, mx] = max(ap);
        [~, mn] = min(ap);

        theseWidths = [theseWidths, abs(mx - mn)/30000]; %#ok<AGROW>
    end

    allWidths = [allWidths, theseWidths]; %#ok<AGROW>
end

if isempty(allWidths)
    error('no spike widths were collected.');
end

%% ---- fit gmm to combined widths ----
numComponents = 2;
gm = fitgmdist(allWidths', numComponents);

means = gm.mu;
stdDevs = sqrt(squeeze(gm.Sigma))';
weights = gm.ComponentProportion;

% sort components so narrow/interneuron first, wide/pyramidal second
[means, order] = sort(means, 'ascend');
stdDevs = stdDevs(order);
weights = weights(order);

intersectionPoint = calculateIntersectionPoint(means, stdDevs, weights);

%% ---- convert to ms for plotting ----
allWidths_ms = allWidths * 1000;
means_ms = means * 1000;
stdDevs_ms = stdDevs * 1000;
intersectionPoint_ms = intersectionPoint * 1000;

%% ---- plot histogram + fitted curves ----
figure('Color','w','Position',[100 100 700 500]);

% slightly smaller bins than before
binWidth_ms = 0.025;

h = histogram(allWidths_ms, 'BinWidth', binWidth_ms, 'EdgeColor', 'black', 'FaceColor', [0.7 0.7 0.7]);
hold on;

x_ms = linspace(min(allWidths_ms), max(allWidths_ms), 1000);

yInt = pdf('Normal', x_ms, means_ms(1), stdDevs_ms(1)) * numel(allWidths_ms) * h.BinWidth * weights(1);
yPyr = pdf('Normal', x_ms, means_ms(2), stdDevs_ms(2)) * numel(allWidths_ms) * h.BinWidth * weights(2);

plot(x_ms, yPyr, 'r', 'LineWidth', 2);
plot(x_ms, yInt, 'g', 'LineWidth', 2);
xline(intersectionPoint_ms, 'k--', 'LineWidth', 2);

xlabel('Spike Width (ms)');
ylabel('Count');
title(sprintf('Combined %s spike widths across %d animals', char(regionToPlot), numFiles));
legend({'Widths','Pyramidal','Interneuron','Intersection'}, 'Location', 'northeast');
box off;
set(gca, 'FontSize', 14, 'LineWidth', 1, 'TickDir', 'out');

fprintf('\ncombined %s results across %d animals:\n', char(regionToPlot), numFiles);
fprintf('n total units = %d\n', numel(allWidths));
fprintf('narrow/interneuron mean = %.6f s (%.3f ms)\n', means(1), means_ms(1));
fprintf('wide/pyramidal mean = %.6f s (%.3f ms)\n', means(2), means_ms(2));
fprintf('intersection point = %.6f s (%.3f ms)\n', intersectionPoint, intersectionPoint_ms);

end

%% ---- helper function for calculating intersection point ----
function intersectionPoint = calculateIntersectionPoint(means, stdDevs, weights)
gaussPDF = @(x,mu,sig) (1/(sig*sqrt(2*pi))) * exp(-(x-mu).^2/(2*sig^2));
intersectionEquation = @(x) ...
    weights(1) * gaussPDF(x,means(1),stdDevs(1)) - ...
    weights(2) * gaussPDF(x,means(2),stdDevs(2));
intersectionPoint = fzero(intersectionEquation, mean(means));
end
