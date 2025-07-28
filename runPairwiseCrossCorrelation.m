function runPairwiseCrossCorrelation(baseDirs)
% pairwise cross-correlation between each interneuron and pyramidal neuron
% code works for both cortex and striatum + saves 3D xcorr matrix (and peak lags and peak correlations)

binSize = 0.001;  % 1 ms bins
maxLagSecs = 1;
maxLagBins = round(maxLagSecs / binSize);
lags = -maxLagBins:maxLagBins;

conslidatedDataFoler = 'X:\David\AnalysesData';
load(fullfile(conslidatedDataFoler, 'AA_classifications.mat'), 'classifications');

% define known base folders (match order used when AA_classifications.mat was created)
animalFolders = {'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData', 'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData', 'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData'};

regions = {'Cortex', 'Striatum'};

allXC = cell(length(baseDirs), length(regions));
allPeakLags = cell(length(baseDirs), length(regions));
allPeakCorrs = cell(length(baseDirs), length(regions));

allXC_Shifted = cell(length(baseDirs), length(regions));
allPeakLags_Shifted = cell(length(baseDirs), length(regions));
allPeakCorrs_Shifted = cell(length(baseDirs), length(regions));

for iRegion = 1:length(regions)
    regionName = regions{iRegion};

    for iDir = 1:length(baseDirs)
        baseDir = baseDirs{iDir};
        fprintf('\nProcessing %s — Session %d: %s\n', regionName, iDir, baseDir);

        neuronFile = fullfile(baseDir, 'neuronDataStruct.mat');
        frFile = fullfile(baseDir, 'NeuralFiringRates10msBins30msGauss.mat');

        if ~isfile(neuronFile) || ~isfile(frFile)
            warning('Missing required files for session %d. Skipping.', iDir);
            continue;
        end

        % Match classification by baseDir
        matchRow = find(contains(animalFolders, baseDir), 1);
        if isempty(matchRow)
            warning('Could not match folderPath to an entry in animalFolders. Skipping.');
            continue;
        end
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

        neuronTypeRow = classifications(matchRow, :);
        regionClass = neuronTypeRow{1, iRegion}(regionInds);
        interFRs = frMatrix(regionClass == 1, :);
        pyrFRs = frMatrix(regionClass == 0, :);

        % shift interneuron matrix uniformly for control
        numBins = size(interFRs, 2);
        minShift = round(30 / binSize);
        maxShift = numBins - minShift;
        rng(42); % or any fixed number as random seed
        shiftAmount = randi([minShift, maxShift]);
        interFRsShifted = circshift(interFRs, [0, shiftAmount]);

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

        fprintf('→ %s — Session %d: %d interneurons × %d pyramidal processed\n', regionName, iDir, numInter, numPyr);

        % compute control %s for each pair (shifted interneurons & unshifted pyramidal)
        numShifts = 100;
        numInter = size(interFRs, 1);
        numPyr = size(pyrFRs, 1);
        numBins = size(interFRs, 2);

        controlCorrs = nan(numInter, numPyr, numShifts);
        minShift = round(30 / binSize);
        maxShift = numBins - minShift;

        for s = 1:numShifts
            shiftAmt = randi([minShift, maxShift]);
            shiftedFRs = circshift(interFRs, [0, shiftAmt]);

            for i = 1:numInter
                for j = 1:numPyr
                    intTS = shiftedFRs(i, :);
                    pyrTS = pyrFRs(j, :);
                    validIdx = ~isnan(intTS) & ~isnan(pyrTS);

                    if sum(validIdx) > 2
                        controlCorrs(i, j, s) = corr(intTS(validIdx)', pyrTS(validIdx)');
                    end
                end
            end
        end

        % calc for percentiles across shifts for each pair
        prc2_5 = prctile(controlCorrs, 2.5, 3);
        prc97_5 = prctile(controlCorrs, 97.5, 3);

    end
end

% save results
save('PairwiseCrossCorrelationResults.mat', 'allXC', 'allPeakLags', 'allPeakCorrs', 'allXC_Shifted', 'allPeakLags_Shifted', 'allPeakCorrs_Shifted', 'lags');

% flatten matrices into vectors
peakLagVec = cell2mat(cellfun(@(x) x(:), allPeakLags(:), 'UniformOutput', false));
peakCorrVec = cell2mat(cellfun(@(x) x(:), allPeakCorrs(:), 'UniformOutput', false));

peakLagVec_shift = cell2mat(cellfun(@(x) x(:), allPeakLags_Shifted(:), 'UniformOutput', false));
peakCorrVec_shift = cell2mat(cellfun(@(x) x(:), allPeakCorrs_Shifted(:), 'UniformOutput', false));

% plot histogram of peak lags
figure;
histogram(peakLagVec, 50, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
hold on;
histogram(peakLagVec_shift, 50, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
xlabel('Peak Lag (sec)');
ylabel('Count');
title('Histogram of Peak Lags');
legend('Unshifted', 'Shifted');
grid on;

% plot histogram of peak correlations
figure;
histogram(peakCorrVec, 50, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
hold on;
histogram(peakCorrVec_shift, 50, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
xlabel('Peak Correlation');
ylabel('Count');
title('Histogram of Peak Correlations');
legend('Unshifted', 'Shifted');
grid on;

end
