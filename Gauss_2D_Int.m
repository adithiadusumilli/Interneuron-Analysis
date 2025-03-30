% Combined Histogram of Widths and Fitted Gaussian Mixture Model for
% D026-032923, D020-062922, D024-111022
% 2D Gaussian with Peak to Valley & Spike Width
% Cortex

neuronDataStructFiles = {'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData\neuronDataStruct.mat', 'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData\neuronDataStruct.mat', 'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData\neuronDataStruct.mat'};
firingRatesFiles = {'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData\NeuralFiringRates10msBins30msGauss.mat', 'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData\NeuralFiringRates10msBins30msGauss.mat', 'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData\NeuralFiringRates10msBins30msGauss.mat'};

allSpikeWidthsCell = cell(1, numel(neuronDataStructFiles));
allCortexWidthsCell = cell(1, numel(neuronDataStructFiles));


for fileIndex = 1:numel(neuronDataStructFiles)
    load(neuronDataStructFiles{fileIndex});
    load(firingRatesFiles{fileIndex});

    cortexWidthsSec = [];
    spikeWidthsSec = [];

    for i = 1:numel(cortexInds)
        waveform = neuronDataStruct(cortexInds(i)).waveforms;
        biggestChan = neuronDataStruct(cortexInds(i)).biggestChan;
        actionPotential = waveform(:, biggestChan);
        [maxVal, maxIndex] = max(actionPotential);
        [minVal, minIndex] = min(actionPotential);
        cortexWidth = abs(maxIndex - minIndex);
        cortexWidthSec = cortexWidth/30000;
        cortexWidthsSec = [cortexWidthsSec, cortexWidthSec];

        if maxVal > abs(minVal)
            peakDirection = 'positive';
        else 
            peakDirection = 'negative';
        end

        if strcmp(peakDirection, 'positive')
            peakHeight = maxVal;
        else
            peakHeight = minVal;
        end
        
        halfPeakHeight = peakHeight/2;
        
        if peakDirection == 'positive'
            leftIndex = find(actionPotential >= halfPeakHeight, 1, 'first');
            rightIndex = find(actionPotential >= halfPeakHeight, 1, 'last');
        else
            leftIndex = find(actionPotential <= halfPeakHeight, 1, 'first');
            rightIndex = find(actionPotential <= halfPeakHeight, 1, 'last');
        end

        spikeWidthSec = abs(rightIndex - leftIndex)/30000;
        spikeWidthsSec = [spikeWidthsSec, spikeWidthSec];
    end

    allCortexWidthsCell{fileIndex} = cortexWidthsSec;
    allSpikeWidthsCell{fileIndex} = spikeWidthsSec;

end

allSpikeWidths = cell2mat(allSpikeWidthsCell);
allCortexWidths = cell2mat(allCortexWidthsCell);

X = [allCortexWidths', allSpikeWidths'];
numComponents = 2;
gm = fitgmdist(X, numComponents);

means = gm.mu;
covariances = gm.Sigma;
proportions = gm.ComponentProportion;

disp('Mean Vectors: ');
disp(means);

figure; 
scatter(allCortexWidths, allSpikeWidths, 'b.');
hold on;

for i = 1:numComponents
    mu = means(i, :);
    sigma = covariances(:, :, i);
    plot_gaussian_contour(mu, sigma);
    %[X1, X2] = meshgrid(linspace(min(allCortexWidths), max(allCortexWidths), 100), linspace(min(allSpikeWidths), max(allSpikeWidths), 100));
   % Z = mvnpdf([X1(:) X2(:)], mu, sigma);
   % Z = reshape(Z, size(X1));
    %contour(X1, X2, Z);
end

xlabel('Cortex Spike Width (seconds)');
ylabel('Half-Peak Width (seconds)');
title('2D Gaussian Mixture Model');
legend('Data Points', 'Gaussian Contours')

xrange = linspace(min(allCortexWidths), max(allCortexWidths), 100);
yrange = linspace(min(allSpikeWidths), max(allSpikeWidths), 100);
[XGrid, YGrid] = meshgrid(xrange, yrange);

Z = zeros(size(XGrid));
for i = 1:numel(xrange)
    for j = 1:numel(yrange)
        xy = [XGrid(j, i), YGrid(j, i)];
        p1 = mvnpdf(xy, means(1, :), covariances(:, :, 1)) * proportions(1);
        p2 = mvnpdf(xy, means(2, :), covariances(:, :, 2)) * proportions(2);
        Z(j, i) = p1 - p2;
    end
end

contour(XGrid, YGrid, Z, [0 0], 'k', 'LineWidth', 2, 'DisplayName', 'Boundary');

annotation('textbox', [0.2, 0.5, 0.1, 0.14], 'String', 'Wide Spiking', 'FitBoxToText', 'on', 'EdgeColor', 'none');
annotation('textbox', [0.6, 0.5, 0.1, 0.1], 'String', 'Narrow Spiking', 'FitBoxToText', 'on', 'EdgeColor', 'none');

%idx_below = Z < 0;
%idx_above = Z > 0;
%mean_coords_below = mean(X(idx_below, :));
%mean_coords_above = mean(X(idx_above, :));
%text(mean_coords_below(1), mean_coords_below(2), 'Wide Spiking', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
%text(mean_coords_below(1), mean_coords_below(2), 'Narrow Spiking', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

hold off;

function plot_gaussian_contour(mu, sigma)
    x = linspace(mu(1) - 3*sqrt(sigma(1,1)), mu(1) + 3*sqrt(sigma(1,1)), 100);
    y = linspace(mu(2) - 3*sqrt(sigma(2,2)), mu(2) + 3*sqrt(sigma(2,2)), 100);
    [X, Y] = meshgrid(x, y);

    Z = zeros(size(X));
    for i = 1:numel(x)
        for j = 1:numel(y)
            Z(j,i) = mvnpdf([X(j, i) Y(j,i)], mu, sigma);
        end
    end

    contour(X, Y, Z, 'Linewidth', 1.5);
end