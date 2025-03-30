% average firing rate barplot
% cortex
% calculates spikes per 10msec (multiply by 100 for spikes/sec)

neuronDataStructFiles = {'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData\neuronDataStruct.mat', 'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData\neuronDataStruct.mat', 'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData\neuronDataStruct.mat'};
firingRatesFiles = {'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData\NeuralFiringRates10msBins30msGauss.mat', 'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData\NeuralFiringRates10msBins30msGauss.mat', 'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData\NeuralFiringRates10msBins30msGauss.mat'};

interneuronFiringRates = [];
pyramidalNeuronFiringRates = [];

for fileIndex = 1:numel(neuronDataStructFiles)
    load(neuronDataStructFiles{fileIndex});
    load(firingRatesFiles{fileIndex});

    for i = 1:numel(cortexInds)
        waveform = neuronDataStruct(cortexInds(i)).waveforms;
        biggestChan = neuronDataStruct(cortexInds(i)).biggestChan;
        actionPotential = waveform(:, biggestChan);
        [maxVal, maxIndex] = max(actionPotential);
        [minVal, minIndex] = min(actionPotential);
        cortexWidth = abs(maxIndex - minIndex);
        cortexWidthsSec = cortexWidth/30000;
        
        timeStamps = neuronDataStruct(cortexInds(i)).timeStamps;
        numSpikes = numel(timeStamps);
        duration = size(cortexFRs, 2);
        firingRate = numSpikes / duration;
       
        if cortexWidthsSec < 0.00036593
            interneuronFiringRates = [interneuronFiringRates, firingRate];
        else
            pyramidalNeuronFiringRates = [pyramidalNeuronFiringRates, firingRate];
        end
    end
end

avgInterneuronFiringRate = mean(interneuronFiringRates);
avgPyramidalNeuronFiringRate = mean(pyramidalNeuronFiringRates);

smeInterneuronFiringRate = std(interneuronFiringRates)/sqrt(length(interneuronFiringRates));
smePyramidalNeuronFiringRate = std(pyramidalNeuronFiringRates)/sqrt(length(pyramidalNeuronFiringRates));

figure;
bar([avgInterneuronFiringRate, avgPyramidalNeuronFiringRate], 'FaceColor', 'flat');
hold on;
errorbar([1, 2], [avgInterneuronFiringRate, avgPyramidalNeuronFiringRate], [smeInterneuronFiringRate, smePyramidalNeuronFiringRate], 'k.', 'LineStyle', 'none', 'CapSize', 17, 'LineWidth', 1.25);
hold off;
xticks([1, 2]);
xticklabels({'Interneurons', 'Pyramidal Neurons'});
ylabel('Average Firing Rate', 'FontSize', 16);
title('Average Firing Rate of Interneurons and Pyramidal Neurons',  'FontSize', 17);
legend('Average Firing Rate', 'SME',  'FontSize', 14);
set(gcf, 'color', 'w');
set(gca, 'LineWidth', 1.5);
set(gca, 'XTicklabel', {'Interneurons', 'Pyramidal Neurons'}, 'FontSize', 15)
box off;

[h, p] = ttest2(interneuronFiringRates, pyramidalNeuronFiringRates);
fprintf('Results of Unpaired T-Test:\n');
fprintf('p-value: %.4f\n', p);

if h 
    fprintf('Statistically significant difference between firing rates.')
else 
    fprintf('No statistically significant difference between firing rates')
end
