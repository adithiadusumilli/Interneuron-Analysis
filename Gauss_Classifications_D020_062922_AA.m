% 1D Gaussian for D020-062922
% Array storing 1s and 0s for pyramidal & interneuron classifications
% Final - Cortex

neuronDataStructFiles = {'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData\neuronDataStruct.mat'};
firingRatesFiles = {'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData\NeuralFiringRates10msBins30msGauss.mat'};

allSpikeWidthsCell = {};

% Load neuron data and compute spike widths
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
        cortexWidthsSec = cortexWidth/30000;  % Convert width to seconds
        cortexWidths = [cortexWidths, cortexWidthsSec];
    end

    allSpikeWidthsCell{fileIndex} = cortexWidths;
end

% Combine all spike widths into one array
allSpikeWidths = cell2mat(allSpikeWidthsCell);

% Fit a Gaussian mixture model
numComponents = 2;
gm = fitgmdist(allSpikeWidths', numComponents);
means = gm.mu;
stdDevs = sqrt(gm.Sigma);
componentProportions = gm.ComponentProportion;

% Plot the histograms and Gaussian curves
figure;
h = histogram(allSpikeWidths, 'BinWidth', 1/20000, 'EdgeColor', 'black', 'FaceColor', 'blue');
hold on;

x = linspace(min(allSpikeWidths), max(allSpikeWidths), 1000);

for i = 1:numComponents
    y = pdf('Normal', x, means(i), stdDevs(i)) * componentProportions(i) * sum(h.Values) * h.BinWidth;
    plot(x, y, 'LineWidth', 2);
end

% Calculate intersection point
intersectionPoint = calculateIntersectionPoint(gm.mu, stdDevs);
disp(['Intersection Point (Width in seconds): ', num2str(intersectionPoint)]);

xlabel('Width (seconds)', 'FontSize', 14);
ylabel('Frequency', 'FontSize', 14);
title('Cortex Spike Widths with Fitted GMMs', 'FontSize', 14);
legend({'Widths', 'Pyramidal', 'Interneurons'});
set(gcf, 'color', 'w');
set(gca, 'LineWidth', 1.5);
box off;
hold off;

% Initialize classification array with NaN for neurons not in cortexInds
classArray = nan(1, numel(neuronDataStruct));

for fileIndex = 1:numel(neuronDataStructFiles)
    load(neuronDataStructFiles{fileIndex});
    load(firingRatesFiles{fileIndex});

    for i = 1:numel(cortexInds)
        cortexWidthSec = allSpikeWidthsCell{fileIndex}(i);
        
        if cortexWidthSec < intersectionPoint
            classArray(cortexInds(i)) = 1; % Interneuron
        else
            classArray(cortexInds(i)) = 0; % Pyramidal neuron
        end
    end

    classifications{fileIndex} = classArray;
end

save('AA_classifications.mat', 'classifications');

disp('Neuron Classifications:');
disp(classifications);

% Print all neurons classified as 1 (interneurons)
for fileIndex = 1:numel(neuronDataStructFiles)
    load(neuronDataStructFiles{fileIndex});
    
    % Get the classification array for this file
    currentClassArray = classifications{fileIndex};
    
    % Find indices of neurons classified as 1 (interneuron)
    interneuronIndices = find(currentClassArray == 1);
    
    % Print the indices of interneurons
    if ~isempty(interneuronIndices)
        fprintf('Interneurons in file %s:\n', neuronDataStructFiles{fileIndex});
        for i = 1:numel(interneuronIndices)
            fprintf('Neuron %d: Classified as Interneuron\n', interneuronIndices(i));
        end
    else
        fprintf('No interneurons found in file %s.\n', neuronDataStructFiles{fileIndex});
    end
    
    % Find indices of neurons classified as 0 (pyramidal neurons)
    pyramidalIndices = find(currentClassArray == 0);
    
    % Print the indices of pyramidal neurons
    if ~isempty(pyramidalIndices)
        fprintf('Pyramidal neurons in file %s:\n', neuronDataStructFiles{fileIndex});
        for i = 1:numel(pyramidalIndices)
            fprintf('Neuron %d: Classified as Pyramidal Neuron\n', pyramidalIndices(i));
        end
    else
        fprintf('No pyramidal neurons found in file %s.\n', neuronDataStructFiles{fileIndex});
    end
end

% Load the neuron classifications
load('AA_classifications.mat'); % Ensure the file path is correct

% Convert classifications from cell array to a numeric array
% Concatenate classifications across all files into a single numeric array
classificationsNumeric = cell2mat(classifications);

% Find the indices for interneurons and pyramidal neurons
interneuronIndices = find(classificationsNumeric == 1);
pyramidalNeuronIndices = find(classificationsNumeric == 0);

% Count the number of interneurons and pyramidal neurons
numInterneurons = numel(interneuronIndices);
numPyramidalNeurons = numel(pyramidalNeuronIndices);

% Calculate the ratio of interneurons to pyramidal neurons
populationRatio = numInterneurons / numPyramidalNeurons;

% Display the counts and ratio
disp(['Number of Interneurons: ', num2str(numInterneurons)]);
disp(['Number of Pyramidal Neurons: ', num2str(numPyramidalNeurons)]);
disp(['Population Ratio (Interneurons to Pyramidal Neurons): ', num2str(populationRatio)]);


% Function to calculate intersection point between two Gaussians
function intersectionPoint = calculateIntersectionPoint(means, stdDevs)
    gaussPDF = @(x, mu, sigma) (1/(sigma * sqrt(2*pi))) * exp(-(x - mu).^2 / (2*sigma^2));
    intersectionEquation = @(x) gaussPDF(x, means(1), stdDevs(1)) - gaussPDF(x, means(2), stdDevs(2));
    intersectionPoint = fzero(intersectionEquation, mean(means));
end
