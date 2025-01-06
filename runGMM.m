function runGMM(baseDirs)
%RUNGMM Summary of this function goes here
%   Detailed explanation goes here

regions = {'Cortex','Striatum'};

for iDir = 1:length(baseDirs)
    allSpikeWidthsCell(iDir,:) = extractSpikeWidths(baseDirs{iDir});
end

for iRegion = 1:length(regions)

    allSpikeWidths = cat(2,allSpikeWidthsCell{:,iRegion});

    numComponents = 2;
    gm = fitgmdist(allSpikeWidths(1,:)', numComponents);
    means = gm.mu;
    stdDevs = sqrt(gm.Sigma);
    colors = {'b', 'r'};
    componentProportions = gm.ComponentProportion;


end

end



function spikeWidths = extractSpikeWidths(baseDir)

neuronDataStructFile = fullfile(baseDir,'neuronDataStruct.mat');
firingRatesFile = fullfile(baseDir,'NeuralFiringRates10msBins30msGauss.mat');

load(neuronDataStructFile)
load(firingRatesFile)

regionsInds = {cortexInds, striatumInds};

for iRegion = 1:length(regionsInds)

    for iNeuron = 1:length(regionsInds{iRegion})

        waveform = neuronDataStruct(regionsInds{iRegion}(iNeuron)).waveforms;
        biggestChan = neuronDataStruct(regionsInds{iRegion}(iNeuron)).biggestChan;
        actionPotential = waveform(:, biggestChan);
        [~, maxIndex] = max(actionPotential);
        [~, minIndex] = min(actionPotential);
        spikeWidth(iNeuron) = abs(maxIndex - minIndex);
        spikeWidthsSec(iNeuron) = spikeWidth(iNeuron)/30000;

    end

    spikeWidths{iRegion} = cat(1,spikeWidth,spikeWidthsSec);

end


end

% 
