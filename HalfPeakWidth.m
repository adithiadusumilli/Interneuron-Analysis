% Combined Histogram of Half Peak Widths in Cortex for D026-032923, D020-062922, D024-111022

neuronDataStructFiles = {'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData\neuronDataStruct.mat', 'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData\neuronDataStruct.mat', 'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData\neuronDataStruct.mat'};
firingRatesFiles = {'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData\NeuralFiringRates10msBins30msGauss.mat', 'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData\NeuralFiringRates10msBins30msGauss.mat', 'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData\NeuralFiringRates10msBins30msGauss.mat'};

allSpikeWidthsCell = {};

for fileIndex = 1:numel(neuronDataStructFiles)
    load(neuronDataStructFiles{fileIndex});
    load(firingRatesFiles{fileIndex});

    spikeWidthsSec = [];

    for i = 1:numel(cortexInds)
        waveform = neuronDataStruct(cortexInds(i)).waveforms;
        biggestChan = neuronDataStruct(cortexInds(i)).biggestChan;
        actionPotential = waveform(:, biggestChan);
        [maxVal, maxIndex] = max(actionPotential);
        [minVal, minIndex] = min(actionPotential);

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
       
        %conditional based on positive or negative so if negative then it
        %should be less than or equal to

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

    allSpikeWidthsCell{fileIndex} = spikeWidthsSec;

end

allSpikeWidths = [];
for fileIndex = 1:numel(allSpikeWidthsCell)
    allSpikeWidths = [allSpikeWidths, allSpikeWidthsCell{fileIndex}];
end

figure;
histogram(allSpikeWidths, 'BinWidth', 1/30000, 'EdgeColor', 'black', 'FaceColor', 'blue');
xlabel('Spike Width at Half Peak (s)');
ylabel('Frequency');
title('Histogram of Half-Peak Widths');