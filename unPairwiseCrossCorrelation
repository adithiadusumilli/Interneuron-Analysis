function runPairwiseCrossCorrelation(baseDirs, classificationFile)
% pairwise cross-correlation between each interneuron and pyramidal neuron
% code works for both cortex and striatum + saves 3D xcorr matrix (and peak lags and peak correlations)

binSize = 0.01;  % 10 ms bins
maxLagSecs = 1;
maxLagBins = round(maxLagSecs / binSize);
lags = -maxLagBins:maxLagBins;

regions = {'Cortex', 'Striatum'};

allXC = cell(length(baseDirs), length(regions));
allPeakLags = cell(length(baseDirs), length(regions));
allPeakCorrs = cell(length(baseDirs), length(regions));

for iRegion = 1:length(regions)
    regionName = regions{iRegion};

    for iDir = 1:length(baseDirs)
        baseDir = baseDirs{iDir};
        fprintf('\nProcessing %s — Session %d: %s\n', regionName, iDir, baseDir);

        classificationPath = fullfile(baseDir, classificationFile);
        neuronFile = fullfile(baseDir, 'neuronDataStruct.mat');
        frFile = fullfile(baseDir, 'NeuralFiringRates10msBins30msGauss.mat');

        if ~isfile(classificationPath) || ~isfile(neuronFile) || ~isfile(frFile)
            warning('Missing required files for session %d. Skipping.', iDir);
            continue;
        end

        load(classificationPath, 'classifications', 'regions');
        load(neuronFile, 'neuronDataStruct');
        load(frFile, 'cortexFRs', 'cortexInds', 'striatumFRs', 'striatumInds');

        % region-specific firing rates and indices
        if iRegion == 1
            frMatrix = cortexFRs;
            regionInds = cortexInds;
        else
            frMatrix = striatumFRs;
            regionInds = striatumInds;
        end

        neuronType = classifications{1, iRegion};
        if isempty(neuronType)
            warning('No classification data for session %d (%s). Skipping.', iDir, regionName);
            continue;
        end

        regionClass = neuronType(regionInds);
        interFRs = frMatrix(regionClass == 1, :);
        pyrFRs = frMatrix(regionClass == 0, :);

        numInter = size(interFRs, 1);
        numPyr = size(pyrFRs, 1);

        if numInter == 0 || numPyr == 0
            warning('No valid neurons in session %d (%s). Skipping.', iDir, regionName);
            continue;
        end

        % setting up for outputs
        xcMat = nan(numInter, numPyr, length(lags));
        peakLagMat = nan(numInter, numPyr);
        peakCorrMat = nan(numInter, numPyr);

        for intIdx = 1:numInter
            intTS = interFRs(intIdx, :);

            for pyrIdx = 1:numPyr
                pyrTS = pyrFRs(pyrIdx, :);

                % manual xcorr & lags
                xc = nan(size(lags));
                for li = 1:length(lags)
                    lag = lags(li);

                    if lag < 0
                        intSeg = intTS(1:end+lag);
                        pyrSeg = pyrTS(1-lag:end);
                    elseif lag > 0
                        intSeg = intTS(1+lag:end);
                        pyrSeg = pyrTS(1:end-lag);
                    else
                        intSeg = intTS;
                        pyrSeg = pyrTS;
                    end

                    % removing nan indices after shifting
                    nanIdx = unique([find(isnan(intSeg)), find(isnan(pyrSeg))]);
                    intSeg(nanIdx) = [];
                    pyrSeg(nanIdx) = [];

                    if length(intSeg) > 2
                        xc(li) = corr(intSeg', pyrSeg');
                    end
                end

                % storing results
                xcMat(intIdx, pyrIdx, :) = xc;
                [peakCorr, peakIdx] = max(xc);
                peakLag = lags(peakIdx) * binSize;

                peakLagMat(intIdx, pyrIdx) = peakLag;
                peakCorrMat(intIdx, pyrIdx) = peakCorr;
            end
        end

        % storing for this session + region
        allXC{iDir, iRegion} = xcMat;
        allPeakLags{iDir, iRegion} = peakLagMat;
        allPeakCorrs{iDir, iRegion} = peakCorrMat;

        fprintf('→ %s — Session %d: %d interneurons × %d pyramidal processed\n', ...
            regionName, iDir, numInter, numPyr);
    end
end

% save results
save('PairwiseCrossCorrelationResults.mat', ...
    'allXC', 'allPeakLags', 'allPeakCorrs', 'lags');

% plot histograms
peakLagVec = cell2mat(cellfun(@(x) x(:), allPeakLags(:), 'UniformOutput', false));
peakCorrVec = cell2mat(cellfun(@(x) x(:), allPeakCorrs(:), 'UniformOutput', false));

figure;
histogram(peakLagVec, 50);
xlabel('Peak Lag (sec)');
ylabel('Count');
title('Histogram of Peak Lags');

figure;
histogram(peakCorrVec, 50);
xlabel('Peak Correlation');
ylabel('Count');
title('Histogram of Peak Correlations');

end
