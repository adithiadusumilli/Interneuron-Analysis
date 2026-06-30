function runCrossCorrelation_saveOnly_getMouseDataNames(mouseIDs, baseSessionNames, probeRegions)
% cortex-only population cross-correlation using getMouseDataNames
% uses new AA_classifications.mat ordering:
% 1 D026, 2 D020, 3 D024, 4 D043, 5 D050, 6 D054

binSize = 0.001;
maxLagSecs = 0.5;

regionName = 'Cortex';
nSess = numel(mouseIDs);

addpath('C:\Users\mirilab\Documents\GlobusTransfer');

if numel(baseSessionNames) ~= nSess || numel(probeRegions) ~= nSess
    error('mouseIDs, baseSessionNames, and probeRegions must all have same length.');
end

peakLags = nan(nSess, 1);
peakCorrs = nan(nSess, 1);
lagCIAll = nan(nSess, 2);

permLagCell = cell(nSess, 1);
animalLabels = mouseIDs(:);

xcorrResults = struct();
xcorrResults.mouseIDs = mouseIDs;
xcorrResults.baseSessionNames = baseSessionNames;
xcorrResults.probeRegions = probeRegions;
xcorrResults.regionName = regionName;
xcorrResults.binSize = binSize;
xcorrResults.maxLagSecs = maxLagSecs;
xcorrResults.sessions = repmat(struct( ...
    'mouseID', '', ...
    'baseSessionName', '', ...
    'processedDataFolder', '', ...
    'animalLabel', '', ...
    'lagsSec', [], ...
    'xc', [], ...
    'peakLag', NaN, ...
    'peakCorr', NaN, ...
    'corrCI', [NaN NaN], ...
    'lagCI', [NaN NaN], ...
    'permPeakLags', []), nSess, 1);

consolidatedDataFolder = 'X:\David\AnalysesData';
load(fullfile(consolidatedDataFolder, 'AA_classifications.mat'), 'classifications');

for iDir = 1:nSess

    mouseID = mouseIDs{iDir};
    baseSessionName = baseSessionNames{iDir};
    probeRegion = probeRegions{iDir};

    fprintf('\nprocessing %s — session %d: %s\n', regionName, iDir, mouseID);

    dataNames = getMouseDataNames(mouseID, baseSessionName, probeRegion);

    frFile = dataNames.NeuralFiringRates1msBins10msGauss;

    if ~isfile(frFile)
        warning('missing firing rate file for %s: %s. skipping.', mouseID, frFile);
        continue;
    end

    load(frFile, 'cortexFRs', 'cortexInds');

    frMatrix = cortexFRs;
    regionInds = cortexInds;

    iRegion = 1; % cortex
    neuronType = classifications{iDir, iRegion};

    if isempty(neuronType)
        warning('no classification data for session %d (%s). skipping.', iDir, mouseID);
        continue;
    end

    regionClass = neuronType(regionInds);

    interneuronFRs = frMatrix(regionClass == 1, :);
    pyramidalFRs   = frMatrix(regionClass == 0, :);

    if isempty(interneuronFRs) || isempty(pyramidalFRs)
        warning('no valid interneuron or pyramidal data in session %d (%s). skipping.', iDir, mouseID);
        continue;
    end

    meanIntRaw = nanmean(interneuronFRs, 1);
    meanPyrRaw = nanmean(pyramidalFRs, 1);

    % real cross-correlation
    [lagsSec, xc, peakLag, peakCorr] = computeManualXCorr(meanIntRaw, meanPyrRaw, binSize, maxLagSecs);

    peakLags(iDir) = peakLag;
    peakCorrs(iDir) = peakCorr;

    % circular-shift control for correlation
    numShifts = 100;
    minShiftBins = round(30 / binSize);
    maxShiftBinsOk = length(meanIntRaw) - minShiftBins;

    controlCorrs = nan(1, numShifts);

    for s = 1:numShifts
        shiftAmt = randi([minShiftBins, maxShiftBinsOk]);
        intShifted = circshift(meanIntRaw, shiftAmt);

        validIdx = ~isnan(intShifted) & ~isnan(meanPyrRaw);
        if sum(validIdx) > 2
            controlCorrs(s) = corr(intShifted(validIdx)', meanPyrRaw(validIdx)');
        end
    end

    prc25 = prctile(controlCorrs, 2.5);
    prc975 = prctile(controlCorrs, 97.5);

    % label permutation null for peak lag
    numPerms = 100;
    permPeakCorrs = nan(1, numPerms);
    permPeakLags = nan(1, numPerms);

    for p = 1:numPerms
        permLabels = regionClass(randperm(numel(regionClass)));

        permIntFRs = frMatrix(permLabels == 1, :);
        permPyrFRs = frMatrix(permLabels == 0, :);

        if isempty(permIntFRs) || isempty(permPyrFRs)
            continue;
        end

        permMeanInt = nanmean(permIntFRs, 1);
        permMeanPyr = nanmean(permPyrFRs, 1);

        [~, ~, permPeakLag, permPeakCorr] = computeManualXCorr(permMeanInt, permMeanPyr, binSize, maxLagSecs);

        permPeakCorrs(p) = permPeakCorr;
        permPeakLags(p) = permPeakLag;
    end

    goodPerms = ~isnan(permPeakCorrs) & ~isnan(permPeakLags);

    if any(goodPerms)
        lagCI = prctile(permPeakLags(goodPerms), [2.5 97.5]);
        permLagCell{iDir} = permPeakLags(goodPerms);
    else
        warning('no valid permutations for session %d; lag CI not computed.', iDir);
        lagCI = [NaN NaN];
        permLagCell{iDir} = [];
    end

    lagCIAll(iDir,:) = lagCI;

    xcorrResults.sessions(iDir).mouseID = mouseID;
    xcorrResults.sessions(iDir).baseSessionName = baseSessionName;
    xcorrResults.sessions(iDir).processedDataFolder = dataNames.processedDataFolder;
    xcorrResults.sessions(iDir).animalLabel = mouseID;
    xcorrResults.sessions(iDir).lagsSec = lagsSec;
    xcorrResults.sessions(iDir).xc = xc;
    xcorrResults.sessions(iDir).peakLag = peakLag;
    xcorrResults.sessions(iDir).peakCorr = peakCorr;
    xcorrResults.sessions(iDir).corrCI = [prc25 prc975];
    xcorrResults.sessions(iDir).lagCI = lagCI;
    xcorrResults.sessions(iDir).permPeakLags = permLagCell{iDir};

    fprintf('→ %s — %s peak lag: %.3f s | peak corr: %.3f\n', ...
        regionName, mouseID, peakLag, peakCorr);
end

fprintf('\n========== summary cortex ==========\n');
fprintf('peak lags (s):\n'); disp(peakLags);
fprintf('peak corrs:\n'); disp(peakCorrs);
fprintf('lag CIs (s):\n'); disp(lagCIAll);

saveDir = fullfile('X:\David\AnalysesData\InterneuronAnalyses\Lab Meeting Pres', ...
    '6 animals run cross correlation, no chunking, pop-wise');

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

savePath = fullfile(saveDir, 'runCrossCorrelation_savedOutputs_all6Animals.mat');

save(savePath, ...
    'xcorrResults', ...
    'peakLags', ...
    'peakCorrs', ...
    'lagCIAll', ...
    'permLagCell', ...
    'animalLabels');

fprintf('\nsaved plotting-ready outputs to:\n%s\n', savePath);

end

function [lagsSec, xc, peakLagSec, peakCorr] = computeManualXCorr(meanIntRaw, meanPyrRaw, binSize, maxLagSecs)

maxLagBins = round(maxLagSecs / binSize);
lagsBins = -maxLagBins:maxLagBins;
xc = nan(size(lagsBins));

for li = 1:length(lagsBins)
    lag = lagsBins(li);

    if lag < 0
        intSeg = meanIntRaw(1:end+lag);
        pyrSeg = meanPyrRaw(1-lag:end);
    elseif lag > 0
        intSeg = meanIntRaw(1+lag:end);
        pyrSeg = meanPyrRaw(1:end-lag);
    else
        intSeg = meanIntRaw;
        pyrSeg = meanPyrRaw;
    end

    validIdx = ~isnan(intSeg) & ~isnan(pyrSeg);

    if sum(validIdx) > 2
        xc(li) = corr(intSeg(validIdx)', pyrSeg(validIdx)');
    end
end

[peakCorr, peakIdx] = max(xc);
lagsSec = lagsBins * binSize;
peakLagSec = lagsSec(peakIdx);

end
