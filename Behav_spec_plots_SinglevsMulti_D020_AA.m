clear;
close all;

% Load Data
load('Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData\neuronDataStruct.mat');
load('Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData\NeuralFiringRates10msBins30msGauss.mat');
load('AA_classifications.mat');

% Scaling
allFRs = allFRs * 100;

% Convert Classifications to Numeric (1/0) if needed
if iscell(classifications)
    classifications = cell2mat(classifications);
end

% Ensure cortexInds are within the valid range (1 to 243)
validNeuronRange = 1:size(allFRs, 1); % 1 to number of neurons (e.g., 1 to 243)
cortexInds = cortexInds(cortexInds <= max(validNeuronRange)); % Remove any indices > 243

% Parameters
numCortexNeurons = numel(cortexInds); % Only use valid cortexInds neurons
numTimepoints = size(allFRs, 2);
numBehaviors = 7; % 7 Behavior Classes
numShifts = 300; % Number of Random Shifts to Create Control Distribution
minShift = 30000;
maxShift = length(regionLabels) - 30000;
zScoreThreshold = 3; % Threshold for Significance (3 Standard Deviations)

% Arrays to Store Modulation Results
modulatedBehaviors = cell(1, numCortexNeurons); % Track which behaviors each neuron is modulated to
multiBehaviorModulated = zeros(1, numCortexNeurons); % Track if neuron is modulated to multiple behaviors

% Counters for Pyramidal and Interneuron Modulation
singleBehaviorPyramidalCount = 0;
multiBehaviorPyramidalCount = 0;
singleBehaviorInterneuronCount = 0;
multiBehaviorInterneuronCount = 0;

% Total Modulated Neurons Counters
totalModulatedPyramidal = 0;
totalModulatedInterneuron = 0;

% Looping through Neurons in cortexInds
for neuronIndex = 1:numCortexNeurons
    actualNeuronIndex = cortexInds(neuronIndex); % Use cortexInds for valid neuron indices
    if actualNeuronIndex > size(allFRs, 1)
        error(['Invalid neuron index: ', num2str(actualNeuronIndex)]);
    end
    
    % Get Firing Rates for Current Neuron
    firingRates = allFRs(actualNeuronIndex, :);
    
    % Calculate Avg Firing Rate per Behavior
    avgFiringRates = NaN(1, numBehaviors); % Initializing with NaNs for Handling Empty Cases
    for behavior = 1:numBehaviors
        % Find Indices for Current Behavior
        behaviorIndices = find(regionLabels == behavior);
        
        % If behaviorIndices is Empty, Skip Calculation
        if isempty(behaviorIndices)
            avgFiringRates(behavior) = NaN; % Set to NaN if no Data for this Behavior
        else
            % Calculate Avg Firing Rate for this Behavior
            behaviorFiringRates = firingRates(behaviorIndices);
            if any(~isnan(behaviorFiringRates))
                avgFiringRates(behavior) = mean(behaviorFiringRates(~isnan(behaviorFiringRates)));
            else
                avgFiringRates(behavior) = NaN; % Ensure NaN if all values are NaN
            end
        end
    end
    
    % Generate Control Distribution by Random Shifts
    shiftedAvgFiringRates = NaN(numShifts, numBehaviors); % Initialize with NaNs to handle empty cases
    for shiftIndex = 1:numShifts
        % Randomly Shift Behavior Labels using Circshift
        shiftAmount = randi([minShift, maxShift]); % Random Shift Amt in ms
        shiftedLabels = circshift(regionLabels, shiftAmount);
        
        % Calculate Avg Firing Rates for Shifted Labels
        for behavior = 1:numBehaviors
            shiftedBehaviorIndices = find(shiftedLabels == behavior);
            
            % If shiftedBehaviorIndices is Empty, Skip Calculation
            if isempty(shiftedBehaviorIndices)
                shiftedAvgFiringRates(shiftIndex, behavior) = NaN; % Set to NaN if no data for this behavior
            else
                % Calculate Avg Firing Rate for Shifted Behavior
                shiftedFiringRates = firingRates(shiftedBehaviorIndices);
                if any(~isnan(shiftedFiringRates))
                    shiftedAvgFiringRates(shiftIndex, behavior) = mean(shiftedFiringRates(~isnan(shiftedFiringRates)));
                else
                    shiftedAvgFiringRates(shiftIndex, behavior) = NaN; % Ensure NaN if all values are NaN
                end
            end
        end
    end
    
    % Compare Unshifted Firing Rates to Control Distribution
    modulatedBehaviors{neuronIndex} = []; % Initialize as empty for current neuron
    for behavior = 1:numBehaviors
        % Extract Control Distribution for Current Behavior
        controlDist = shiftedAvgFiringRates(:, behavior);
        
        % Check for NaNs or Zero STDev
        if all(isnan(controlDist)) || std(controlDist) == 0
            continue;
        end
        
        % Calculate Z-score of Original Firing Rate Compared to Control Distribution
        originalFiringRate = avgFiringRates(behavior);
        if isnan(originalFiringRate)
            continue;
        end
        zScore = (originalFiringRate - mean(controlDist)) / std(controlDist);
        
        % Determine Modulation Based on Z-Score Threshold
        if zScore >= zScoreThreshold
            % Track the modulated behavior for the neuron
            modulatedBehaviors{neuronIndex} = [modulatedBehaviors{neuronIndex}, behavior];
        end
    end
    
    % If neuron is modulated, classify it based on behavior count
    if ~isempty(modulatedBehaviors{neuronIndex})
        if numel(modulatedBehaviors{neuronIndex}) > 1
            multiBehaviorModulated(neuronIndex) = 1; % Neuron is modulated to multiple behaviors
        end
        
        % Check if Neuron is Pyramidal or Interneuron
        neuronType = '';
        if classifications(actualNeuronIndex) == 0
            neuronType = 'Pyramidal';
            totalModulatedPyramidal = totalModulatedPyramidal + 1;
            if multiBehaviorModulated(neuronIndex)
                multiBehaviorPyramidalCount = multiBehaviorPyramidalCount + 1;
            else
                singleBehaviorPyramidalCount = singleBehaviorPyramidalCount + 1;
            end
        else
            neuronType = 'Interneuron';
            totalModulatedInterneuron = totalModulatedInterneuron + 1;
            if multiBehaviorModulated(neuronIndex)
                multiBehaviorInterneuronCount = multiBehaviorInterneuronCount + 1;
            else
                singleBehaviorInterneuronCount = singleBehaviorInterneuronCount + 1;
            end
        end
        
        % Plot Avg Firing Rates per Behavior for the Current Neuron
        figure;
        set(gcf, 'Color', 'w'); % Set figure background to white
        set(gca, 'Color', 'w'); % Set axes background to white
        plot(shiftedAvgFiringRates', 'color', [0.7, 0.7, 0.7], 'LineWidth', 0.5); % Shifted
        hold on;
        plot(1:numBehaviors, avgFiringRates, '-o', 'MarkerFaceColor', 'auto'); % Unshifted, filled circles
        
        title(['Neuron ', num2str(neuronIndex), ' (Cortex Index ', num2str(cortexInds(neuronIndex)), ') - ', neuronType]);
        xlabel('Behavior');
        ylabel('Average Firing Rate (Hz)');
        xticks(1:numBehaviors);
        xticklabels({'Climb Up', 'Climb Down', 'Misc/Jump', 'Walk', 'Misc/Rear/Still', 'Groom', 'Eat'});
        grid on;

        % Display modulated behaviors for the neuron
        disp(['Neuron ', num2str(neuronIndex), ' (Cortex Index ', num2str(cortexInds(neuronIndex)), ') is modulated to behaviors: ', num2str(modulatedBehaviors{neuronIndex})]);
        if multiBehaviorModulated(neuronIndex)
            disp(['Neuron ', num2str(neuronIndex), ' is modulated to multiple behaviors']);
        else
            disp(['Neuron ', num2str(neuronIndex), ' is modulated to a single behavior']);
        end
    end
end

% Calculate Total Number of Modulated Neurons
totalModulatedPyramidal = singleBehaviorPyramidalCount + multiBehaviorPyramidalCount;
totalModulatedInterneuron = singleBehaviorInterneuronCount + multiBehaviorInterneuronCount;

% Calculate Percentages
if totalModulatedPyramidal > 0
    pyramidalSinglePercentage = (singleBehaviorPyramidalCount / totalModulatedPyramidal) * 100;
    pyramidalMultiPercentage = (multiBehaviorPyramidalCount / totalModulatedPyramidal) * 100;
else
    pyramidalSinglePercentage = 0;
    pyramidalMultiPercentage = 0;
end

if totalModulatedInterneuron > 0
    interneuronSinglePercentage = (singleBehaviorInterneuronCount / totalModulatedInterneuron) * 100;
    interneuronMultiPercentage = (multiBehaviorInterneuronCount / totalModulatedInterneuron) * 100;
else
    interneuronSinglePercentage = 0;
    interneuronMultiPercentage = 0;
end

% Data for Bar Plot
data = [
    pyramidalSinglePercentage, pyramidalMultiPercentage;
    interneuronSinglePercentage, interneuronMultiPercentage
];

% Save modulation counts and percentages to a .mat file
modulationCounts.singleBehaviorPyramidalCount = singleBehaviorPyramidalCount;
modulationCounts.multiBehaviorPyramidalCount = multiBehaviorPyramidalCount;
modulationCounts.singleBehaviorInterneuronCount = singleBehaviorInterneuronCount;
modulationCounts.multiBehaviorInterneuronCount = multiBehaviorInterneuronCount;
modulationCounts.totalModulatedPyramidal = totalModulatedPyramidal;
modulationCounts.totalModulatedInterneuron = totalModulatedInterneuron;
modulationCounts.pyramidalSinglePercentage = pyramidalSinglePercentage;
modulationCounts.pyramidalMultiPercentage = pyramidalMultiPercentage;
modulationCounts.interneuronSinglePercentage = interneuronSinglePercentage;
modulationCounts.interneuronMultiPercentage = interneuronMultiPercentage;

% Save the counts to a .mat file
save('modulationCounts.mat', 'modulationCounts');

% Create Bar Plot
figure;
set(gcf, 'Color', 'w'); % Set figure background to white
set(gca, 'Color', 'w'); % Set axes background to white
bar(data, 'grouped');
set(gca, 'XTickLabel', {'Pyramidal Neurons', 'Interneurons'});
ylabel('Percentage of Modulated Neurons');
legend('Single Behavior Modulated', 'Multi Behavior Modulated');
title('Percentage of Neurons Modulated to Single vs Multiple Behaviors');
grid on;

