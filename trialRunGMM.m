function trialRunGMM(baseDirs)
% function runs GMM and classifies neurons for Cortex and Striatum

regions = {'Cortex', 'Striatum'};
numRegions = length(regions);

allSpikeWidthsCell = cell(length(baseDirs), numRegions);
allNeuronIndices = cell(length(baseDirs), numRegions);
classifications = cell(length(baseDirs), numRegions);

% extract spike widths and neuron indices from each session
for iDir = 1:length(baseDirs)
    [spikeWidthsPerRegion, neuronIndsPerRegion] = extractSpikeWidths(baseDirs{iDir});
    for iRegion = 1:numRegions
        allSpikeWidthsCell{iDir, iRegion} = spikeWidthsPerRegion{iRegion};
        allNeuronIndices{iDir, iRegion} = neuronIndsPerRegion{iRegion};
    end
end

% GMM fitting and classification
for iRegion = 1:numRegions
    regionName = regions{iRegion};
    allWidths = cell2mat(allSpikeWidthsCell(:, iRegion)');

    % skip if no data
    if isempty(allWidths)
        warning('%s has no spike width data. Skipping...', regionName);
        continue;
    end

    % fit 2-component GMM
    gm = fitgmdist(allWidths', 2);
    means = gm.mu;
    stdDevs = sqrt(squeeze(gm.Sigma));
    proportions = gm.ComponentProportion;

    % plot GMM
    % first sort components by mean to assign labels
    [sortedMeans, sortIdx] = sort(means);
    sortedStdDevs = stdDevs(sortIdx);
    sortedProportions = proportions(sortIdx);
    labels = {'Interneurons', 'Pyramidal Neurons'};

    % plot fig
    figure('Name', [regionName ' Spike Widths'], 'NumberTitle', 'off', 'Visible', 'on');
    h = histogram(allWidths, 'BinWidth', 1/20, 'EdgeColor', 'black', 'FaceColor', [0.7 0.7 0.9], 'DisplayName', '', 'HandleVisibility', 'off');
    hold on;

    x = linspace(min(allWidths), max(allWidths), 1000);
    colors = lines(3);
    for region = 1:2
        y = pdf('Normal', x, sortedMeans(region), sortedStdDevs(region)) * ...
            sortedProportions(region) * sum(h.Values) * h.BinWidth;
        plot(x, y, 'LineWidth', 2, 'DisplayName', labels{region},'Color',colors(region,:));
    end

    xlabel('Spike Width (seconds)', 'FontSize', 14);
    ylabel('Frequency', 'FontSize', 14);
    title([regionName ' Spike Widths with GMM'], 'FontSize', 14);
    legend show;
    set(gcf, 'Color', 'w');
    box off;
    hold off;
    drawnow;

    % Calculate and report intersection point
    intersection = calculateIntersectionPoint(means, stdDevs);
    fprintf('%s Intersection Point: %.6f seconds\n', regionName, intersection);

    % Classify neurons
    for iDir = 1:length(baseDirs)
        spikeWidths = allSpikeWidthsCell{iDir, iRegion};
        neuronInds = allNeuronIndices{iDir, iRegion};
        neuronType = nan(1, max(neuronInds));

        for i = 1:length(spikeWidths)
            if spikeWidths(i) < intersection
                neuronType(neuronInds(i)) = 1;  % Interneuron
            else
                neuronType(neuronInds(i)) = 0;  % Pyramidal
            end
        end

        classifications{iDir, iRegion} = neuronType;
    end
end

% Save results
save('GMM_classifications.mat', 'classifications', 'regions', 'baseDirs');

end

function [spikeWidths, neuronIndices] = extractSpikeWidths(baseDir)
% function extracts spike widths and neuron indices for cortex and striatum

neuronDataStructFile = fullfile(baseDir, 'neuronDataStruct.mat');
firingRatesFile = fullfile(baseDir, 'NeuralFiringRates10msBins30msGauss.mat');

if ~isfile(neuronDataStructFile)
    error('File not found: %s', neuronDataStructFile);
end
if ~isfile(firingRatesFile)
    error('File not found: %s', firingRatesFile);
end

load(neuronDataStructFile);
load(firingRatesFile);

regionsInds = {cortexInds, striatumInds};
spikeWidths = cell(1, 2);
neuronIndices = cell(1, 2);

for iRegion = 1:2
    inds = regionsInds{iRegion};
    widthsSec = nan(1, numel(inds));

    for i = 1:numel(inds)
        neuron = neuronDataStruct(inds(i));
        waveform = neuron.waveforms(:, neuron.biggestChan);
        [~, maxIdx] = max(waveform);
        [~, minIdx] = min(waveform);
        widthsSec(i) = abs(maxIdx - minIdx) / 30;
    end

    spikeWidths{iRegion} = widthsSec;
    neuronIndices{iRegion} = inds;
end

end

function intersectionPoint = calculateIntersectionPoint(means, stdDevs)
% function (same as original) calculates intersection between two Gaussians

gaussPDF = @(x, mu, sigma) (1 / (sigma * sqrt(2*pi))) * ...
    exp(-(x - mu).^2 / (2 * sigma^2));

intersectionEquation = @(x) gaussPDF(x, means(1), stdDevs(1)) - ...
                            gaussPDF(x, means(2), stdDevs(2));

intersectionPoint = fzero(intersectionEquation, mean(means));
end