function runCrossCorrelation_saveOnly_getMouseDataNames(mouseIDs, baseSessionNames, probeRegions)
% cortex-only population cross-correlation using getMouseDataNames
% uses new AA_classifications.mat ordering:
% 1 D026, 2 D020, 3 D024, 4 D043, 5 D050, 6 D054

% updated:
% label permutation null now uses the same perm-threshold idea as chunked:
%   shiftCorrUpper95 = 95th percentile of circular-shift control correlations
%   keep drawing label permutations until permPeakCorr >= shiftCorrUpper95
%   or until maxPermTries is reached

% 300 accepted-permutation slots are generated:
%   permutations   1:100 = raw H0 source lags
%   permutations 101:200 = raw H+50 source lags
%   permutations 201:300 = raw H-50 source lags

% All 300 peak lags are saved unshifted. The later post hoc analysis applies
% the model-specific +50 ms and -50 ms transformations.

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

numShifts = 100;
numPermsPerHypothesis = 100;
numPerms = 3 * numPermsPerHypothesis;  % 300 total permutation slots
maxPermTries = 1000;

permIndsH0 = 1:numPermsPerHypothesis;
permIndsH50 = (numPermsPerHypothesis + 1):(2 * numPermsPerHypothesis);
permIndsHneg50 = (2 * numPermsPerHypothesis + 1):(3 * numPermsPerHypothesis);

%% ---------------- initialize results structure ----------------

xcorrResults = struct();

xcorrResults.mouseIDs = mouseIDs;
xcorrResults.baseSessionNames = baseSessionNames;
xcorrResults.probeRegions = probeRegions;
xcorrResults.regionName = regionName;
xcorrResults.binSize = binSize;
xcorrResults.maxLagSecs = maxLagSecs;

xcorrResults.numShifts = numShifts;
xcorrResults.numPerms = numPerms;
xcorrResults.maxPermTries = maxPermTries;

xcorrResults.numPermsPerHypothesis = numPermsPerHypothesis;
xcorrResults.permIndsH0 = permIndsH0;
xcorrResults.permIndsH50 = permIndsH50;
xcorrResults.permIndsHneg50 = permIndsHneg50;

xcorrResults.sessions = repmat(struct('mouseID', '', 'baseSessionName', '', 'processedDataFolder', '', 'animalLabel', '', 'lagsSec', [], ...
    'xc', [], 'peakLag', NaN, 'peakCorr', NaN, 'corrCI', [NaN NaN], 'shiftCorrUpper95', NaN, 'controlCorrs', [], 'lagCI', [NaN NaN], ...
    'permPeakLags', [], 'permPeakLagsAll300', [], 'permPeakCorrs', [], 'permTryCounts', [], 'permAccepted', [], 'permIndsH0', [], ...
    'permIndsH50', [], 'permIndsHneg50', []), nSess, 1);

%% ---------------- load classifications ----------------

consolidatedDataFolder = 'X:\David\AnalysesData';

load( ...
    fullfile(consolidatedDataFolder, 'AA_classifications.mat'), ...
    'classifications');

%% ========================================================================
%% loop through sessions
%% ========================================================================

for iDir = 1:nSess

    mouseID = mouseIDs{iDir};
    baseSessionName = baseSessionNames{iDir};
    probeRegion = probeRegions{iDir};

    fprintf('\nprocessing %s — session %d: %s\n', regionName, iDir, mouseID);

    dataNames = getMouseDataNames( mouseID, baseSessionName, probeRegion);

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
    pyramidalFRs = frMatrix(regionClass == 0, :);

    if isempty(interneuronFRs) || isempty(pyramidalFRs)
        warning('no valid interneuron or pyramidal data in session %d (%s). skipping.', iDir, mouseID);
        continue;
    end

    meanIntRaw = nanmean(interneuronFRs, 1);
    meanPyrRaw = nanmean(pyramidalFRs, 1);

    %% ---------------- real cross-correlation ----------------

    [lagsSec, xc, peakLag, peakCorr] = computeManualXCorr(meanIntRaw, ...
        meanPyrRaw, binSize, maxLagSecs);

    peakLags(iDir) = peakLag;
    peakCorrs(iDir) = peakCorr;

    %% ---------------- circular-shift correlation control ----------------

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

    goodShift = ~isnan(controlCorrs) & isfinite(controlCorrs);
    if any(goodShift)
        prc25 = prctile( controlCorrs(goodShift), 2.5);
        prc975 = prctile(controlCorrs(goodShift), 97.5);
        shiftCorrUpper95 = prctile(controlCorrs(goodShift),95);
    else

        prc25 = NaN;
        prc975 = NaN;
        shiftCorrUpper95 = NaN;
    end

    fprintf('shiftCorrUpper95 for %s = %.6f\n', mouseID, shiftCorrUpper95);

    %% ---------------- label permutation null ----------------
    % Label permutation null for peak lag with permutation-threshold filtering.
    
    % The arrays below remain ordered and length 300 so that:
    %   1:100   = H0 source permutations
    %   101:200 = H+50 source permutations
    %   201:300 = H-50 source permutations

    permPeakCorrs = nan(1, numPerms);
    permPeakLags = nan(1, numPerms);
    permTryCounts = nan(1, numPerms);
    permAccepted = false(1, numPerms);

    for p = 1:numPerms
        tryCount = 0;
        acceptedPerm = false;

        while ~acceptedPerm && tryCount < maxPermTries
            tryCount = tryCount + 1;
            permLabels = regionClass( ...
                randperm(numel(regionClass)));
            permIntFRs = frMatrix(permLabels == 1, :);
            permPyrFRs = frMatrix(permLabels == 0, :);

            if isempty(permIntFRs) || isempty(permPyrFRs)
                continue;
            end

            permMeanInt = nanmean(permIntFRs, 1);
            permMeanPyr = nanmean(permPyrFRs, 1);

            [~, ~, permPeakLag, permPeakCorr] = computeManualXCorr(permMeanInt, permMeanPyr, binSize, maxLagSecs);

            % Match chunked code logic: accept permutation if its peak correlation reaches the 95th-percentile circular-shift correlation threshold.
            if isnan(shiftCorrUpper95) || permPeakCorr >= shiftCorrUpper95
                acceptedPerm = true;
                permPeakCorrs(p) = permPeakCorr;
                permPeakLags(p) = permPeakLag;
            end
        end

        permTryCounts(p) = tryCount;
        permAccepted(p) = acceptedPerm;

        if ~acceptedPerm
            warning('Permutation %d for %s was not accepted after %d tries.', p, mouseID, maxPermTries);
        end
    end

    %% ---------------- define valid permutations by hypothesis block ----------------

    goodPermsAll = ~isnan(permPeakCorrs) & ~isnan(permPeakLags) & permAccepted;

    goodPermsH0 = false(1, numPerms);
    goodPermsH0(permIndsH0) = goodPermsAll(permIndsH0);

    goodPermsH50 = false(1, numPerms);
    goodPermsH50(permIndsH50) = goodPermsAll(permIndsH50);

    goodPermsHneg50 = false(1, numPerms);
    goodPermsHneg50(permIndsHneg50) = ...
        goodPermsAll(permIndsHneg50);

    %% ---------------- H0 permutation lag CI ----------------
    % Keep the existing permutation lag CI based only on the first
    % 100 H0 permutations.

    if any(goodPermsH0)
        lagCI = prctile(permPeakLags(goodPermsH0),[2.5 97.5]);
        permLagCell{iDir} = permPeakLags(goodPermsH0);
    else

        warning(['no valid accepted H0 permutations for session %d; ' 'lag CI not computed.'], iDir);
        lagCI = [NaN NaN];
        permLagCell{iDir} = [];
    end

    lagCIAll(iDir,:) = lagCI;

    %% ---------------- store results ----------------

    xcorrResults.sessions(iDir).mouseID = mouseID;
    xcorrResults.sessions(iDir).baseSessionName = baseSessionName;
    xcorrResults.sessions(iDir).processedDataFolder = dataNames.processedDataFolder;
    xcorrResults.sessions(iDir).animalLabel = mouseID;
    xcorrResults.sessions(iDir).lagsSec = lagsSec;
    xcorrResults.sessions(iDir).xc = xc;
    xcorrResults.sessions(iDir).peakLag = peakLag;
    xcorrResults.sessions(iDir).peakCorr = peakCorr;
    xcorrResults.sessions(iDir).corrCI = [prc25 prc975];
    xcorrResults.sessions(iDir).shiftCorrUpper95 = shiftCorrUpper95;
    xcorrResults.sessions(iDir).controlCorrs = controlCorrs;
    xcorrResults.sessions(iDir).lagCI = lagCI;
    
    % Legacy H0-only valid lag vector retained for the current plotting functions.
    xcorrResults.sessions(iDir).permPeakLags = permLagCell{iDir};

    % Ordered full 300-element lag array for the post hoc model analysis. NaNs are retained so the three index blocks remain fixed.
    xcorrResults.sessions(iDir).permPeakLagsAll300 = permPeakLags;

    % These are already full ordered 300-element arrays, so no duplicate
    % All300 copies are necessary
   
    xcorrResults.sessions(iDir).permPeakCorrs = permPeakCorrs;
    xcorrResults.sessions(iDir).permTryCounts = permTryCounts;
    xcorrResults.sessions(iDir).permAccepted = permAccepted;
    xcorrResults.sessions(iDir).permIndsH0 = permIndsH0;
    xcorrResults.sessions(iDir).permIndsH50 = permIndsH50;
    xcorrResults.sessions(iDir).permIndsHneg50 = permIndsHneg50;

    %% ---------------- print session results ----------------

    fprintf(['→ %s — %s peak lag: %.3f s | peak corr: %.3f | ' 'accepted perms: %d/%d | median tries: %.1f\n'], regionName, mouseID, peakLag, peakCorr, ...
        sum(permAccepted), numPerms, median(permTryCounts, 'omitnan'));

    fprintf(['  accepted by block: H0=%d/%d | H+50=%d/%d | ' 'H-50=%d/%d\n'], sum(permAccepted(permIndsH0)), ...
        numPermsPerHypothesis, sum(permAccepted(permIndsH50)), numPermsPerHypothesis, sum(permAccepted(permIndsHneg50)), ...
        numPermsPerHypothesis);
end

%% ---------------- summary ----------------

fprintf('\n========== summary cortex ==========\n');

fprintf('peak lags (s):\n');
disp(peakLags);

fprintf('peak corrs:\n');
disp(peakCorrs);

fprintf('lag CIs (s):\n');
disp(lagCIAll);

%% ---------------- save ----------------

saveDir = fullfile('X:\David\AnalysesData\InterneuronAnalyses\Lab Meeting Pres', '6 animals run cross correlation, no chunking, pop-wise');

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

savePath = fullfile(saveDir, 'runCrossCorrelation_savedOutputs_all6Animals.mat');
save(savePath, 'xcorrResults', 'peakLags', 'peakCorrs', 'lagCIAll', 'permLagCell', 'animalLabels', '-v7.3');
fprintf('\nsaved plotting-ready outputs to:\n%s\n', savePath);

end

%% ========================================================================
%% helper: manual cross-correlation
%% ========================================================================

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
