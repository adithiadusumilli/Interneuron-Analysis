% Combined Histogram of Widths and Fitted Gaussian Mixture Model for
% D026-032923, D020-062922, D024-111022
% Final - Cortex

neuronDataStructFiles = {'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData\neuronDataStruct.mat', 'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData\neuronDataStruct.mat', 'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData\neuronDataStruct.mat'};
firingRatesFiles = {'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData\NeuralFiringRates10msBins30msGauss.mat', 'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData\NeuralFiringRates10msBins30msGauss.mat', 'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData\NeuralFiringRates10msBins30msGauss.mat'};

allSpikeWidthsCell = {};

for fileIndex = 1:numel(neuronDataStructFiles)
    load(neuronDataStructFiles{fileIndex});
    load(firingRatesFiles{fileIndex});

    cortexWidths = [];

    for i = 1:numel(cortexInds)
        waveform = neuronDataStruct(cortexInds(i)).waveforms;
        biggestChan = neuronDataStruct(cortexInds(i)).biggestChan;
        actionPotential = waveform(:, biggestChan);
        [maxVal, maxIndex] = max(actionPotential);
        [minVal, minIndex] = min(actionPotential);
        cortexWidth = abs(maxIndex - minIndex);
        cortexWidthsSec = cortexWidth/30000;
        cortexWidths = [cortexWidths, cortexWidthsSec];
    end

    allSpikeWidthsCell{fileIndex} = cortexWidths;

end

allSpikeWidths = cell2mat(allSpikeWidthsCell);


numComponents = 2;
gm = fitgmdist(allSpikeWidths', numComponents);
means = gm.mu;
stdDevs = sqrt(gm.Sigma);
colors = {'b', 'r'};
componentProportions = gm.ComponentProportion;

[counts, binCenters] = hist(allSpikeWidths, 50);
[~, peaks] = findpeaks(counts, 'SortStr', 'descend', 'NPeaks', numComponents);
%peaksIdxSorted = sort(idx);

figure;
h = histogram(allSpikeWidths, 'BinWidth', 1/20000, 'EdgeColor', 'black', 'FaceColor', 'blue');
hold on;

x = linspace(min(allSpikeWidths), max(allSpikeWidths), 1000);

for i = 1:numComponents
    y = pdf('Normal', x, means(i), stdDevs(i)) * componentProportions(i) * sum(h.Values) * h.BinWidth;
    plot(x, y, 'Color', colors{i}, 'LineWidth', 2);
    %pdf_i = pdf('Normal', x, means(i), stdDevs(i));
    %scaleFactor = peaks(i)/max(pdf_i);
    %pdf_i = pdf_i * scaleFactor;
    % y = y + pdf_i * sum(h.Values) * h.BinWidth;
end

intersectionPoint = calculateIntersectionPoint(gm.mu, stdDevs);
disp(['Intersection Point (Width in seconds): ', num2str(intersectionPoint)]);

% Plot the histogram of all spike widths with blue color and thinner bars
figure;
h = histogram(allSpikeWidths, 'BinWidth', 1/20000, 'EdgeColor', 'black', 'FaceColor', 'blue'); % Set FaceColor to blue and adjust BinWidth for thinner bars
hold on;

% Define x-axis range for Gaussian curve plots
x = linspace(min(allSpikeWidths), max(allSpikeWidths), 1000);

% Plot Gaussian Mixture Model curves for pyramidal and interneuron distributions
pyramidal_curve = plot(x, pdf('Normal', x, means(1), stdDevs(1)) * componentProportions(1) * sum(h.Values) * h.BinWidth, ...
                       'b-', 'LineWidth', 2); % Blue curve for pyramidal neurons
interneuron_curve = plot(x, pdf('Normal', x, means(2), stdDevs(2)) * componentProportions(2) * sum(h.Values) * h.BinWidth, ...
                         'r-', 'LineWidth', 2); % Red curve for interneurons

% Display intersection point
intersectionPoint = calculateIntersectionPoint(gm.mu, stdDevs);
disp(['Intersection Point (Width in seconds): ', num2str(intersectionPoint)]);

% Set labels and title
xlabel('Voltage Durations (s)', 'FontSize', 20);
ylabel('Number of Neurons', 'FontSize', 20);
title('M1 Voltage Durations with Fit Gaussian Mixture Models', 'FontSize', 22);

% Only include Gaussian curves in the legend, excluding the histogram bars
legend([pyramidal_curve, interneuron_curve], {'Pyramidal Fit', 'Interneuron Fit'}, 'FontSize', 14, 'Location', 'best');

% Set figure and axis properties
set(gcf, 'color', 'w');
set(gca, 'LineWidth', 1.5, 'FontSize', 16); % Increase axis tick label size

% Final adjustments
box off; % Remove top and right axis lines
hold off;

function intersectionPoint = calculateIntersectionPoint(means, stdDevs)
    gaussPDF = @(x, mu, sigma) (1/(sigma * sqrt(2*pi))) * exp(-(x - mu).^2 / (2*sigma^2));
    intersectionEquation = @(x) gaussPDF(x, means(1), stdDevs(1)) - gaussPDF(x, means(2), stdDevs(2));
    intersectionPoint = fzero(intersectionEquation, mean(means));
end

    