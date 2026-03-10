function runCrossCorrelation_saveOnly(baseDirs)
% runCrossCorrelation_saveOnly – cortex-only cross-correlation + permutation lag ci

% - cortex only (uses cortexFRs / cortexInds and iRegion = 1 in AA_classifications)
% - computes interneuron vs pyramidal population xc for each session
% - builds null dist of:
%     (a) corr via circular shifts
%     (b) peak lag via label perms
% - saves plotting-ready outputs so figures can be regenerated without rerunning xcorr

% j run:
% baseDirs = {
%     'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData', ...
%     'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData', ...
%     'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData', ...
%     'X:\David\ArenaRecordings\D043-020525-ArenaRecording\ProcessedData'
% };
% runCrossCorrelation_saveOnly(baseDirs)

binSize = 0.001;   % 1 ms bins
maxLagSecs = 0.5;
maxLagBins = round(maxLagSecs / binSize); %#ok<NASGU>

regionName = 'Cortex';
nSess = numel(baseDirs);

peakLags = nan(nSess, 1);
peakCorrs = nan(nSess, 1);
lagCIAll = nan(nSess, 2);

permLagCell = cell(nSess, 1);
animalLabels = cell(nSess, 1);

% store all plotting-ready outputs
xcorrResults = struct();
xcorrResults.baseDirs = baseDirs;
xcorrResults.regionName = regionName;
xcorrResults.binSize = binSize;
xcorrResults.maxLagSecs = maxLagSecs;
xcorrResults.sessions = repmat(struct( ...
    'baseDir', '', ...
    'animalLabel', '', ...
    'lagsSec', [], ...
    'xc', [], ...
    'peakLag', NaN, ...
    'peakCorr', NaN, ...
    'corrCI', [NaN NaN], ...
    'lagCI', [NaN NaN], ...
    'permPeakLags', []), nSess, 1);

conslidatedDataFoler = 'X:\David\AnalysesData';
load(fullfile(conslidatedDataFoler, 'AA_classifications.mat'), 'classifications');

animalFolders = {
    'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData', ...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData', ...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData', ...
    'X:\David\ArenaRecordings\D043-020425-ArenaRecording\ProcessedData'
};

for iDir = 1:nSess
    baseDir = baseDirs{iDir};
    fprintf('\nprocessing %s — session %d: %s\n', regionName, iDir, baseDir);

    [~, folderName] = fileparts(baseDir);
    tag = regexp(folderName, 'D\d+', 'match', 'once');
    if isempty(tag)
        animalLabels{iDir} = folderName;
    else
        animalLabels{iDir} = tag;
    end

    neuronFile = fullfile(baseDir, 'neuronDataStruct.mat');
    frFile = fullfile(baseDir, 'NeuralFiringRates1msBins10msGauss.mat');

    if ~isfile(neuronFile) || ~isfile(frFile)
        warning('missing files in %s. skipping.', baseDir);
        continue;
    end

    load(neuronFile, 'neuronDataStruct'); %#ok<NASGU>
    load(frFile, 'cortexFRs', 'cortexInds');

    frMatrix = cortexFRs;
    regionInds = cortexInds;

    matchRow = find(contains(animalFolders, baseDir), 1);
    if isempty(matchRow)
        warning('could not match baseDir to animalFolders list. skipping.');
        continue;
    end

    iRegion = 1; % cortex
    neuronType = classifications{matchRow, iRegion};

    if isempty(neuronType)
        warning('no classification data for session %d (%s). skipping.', iDir, regionName);
        continue;
    end

    regionClass = neuronType(regionInds);
    interneuronFRs = frMatrix(regionClass == 1, :);
    pyramidalFRs = frMatrix(regionClass == 0, :);

    if isempty(interneuronFRs) || isempty(pyramidalFRs)
        warning('no valid interneuron or pyramidal data in session %d (%s). skipping.', iDir, regionName);
        continue;
    end

    meanIntRaw = nanmean(interneuronFRs, 1);
    meanPyrRaw = nanmean(pyramidalFRs, 1);

    % ===== non-permuted xc =====
    [lagsSec, xc, peakLag, peakCorr] = computeManualXCorr(meanIntRaw, meanPyrRaw, binSize, maxLagSecs);
    peakLags(iDir)  = peakLag;
    peakCorrs(iDir) = peakCorr;

    % ===== control (circular shift) distribution for correlation =====
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

    % ===== permutation test for lag =====
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
        warning('no valid permutations for session %d; lag ci not computed.', iDir);
        lagCI = [NaN NaN];
        permLagCell{iDir} = [];
    end

    lagCIAll(iDir, :) = lagCI;

    xcorrResults.sessions(iDir).baseDir = baseDir;
    xcorrResults.sessions(iDir).animalLabel = animalLabels{iDir};
    xcorrResults.sessions(iDir).lagsSec = lagsSec;
    xcorrResults.sessions(iDir).xc = xc;
    xcorrResults.sessions(iDir).peakLag = peakLag;
    xcorrResults.sessions(iDir).peakCorr = peakCorr;
    xcorrResults.sessions(iDir).corrCI = [prc25 prc975];
    xcorrResults.sessions(iDir).lagCI = lagCI;
    xcorrResults.sessions(iDir).permPeakLags = permLagCell{iDir};

    fprintf('→ %s — session %d  peak lag: %.3f s | peak corr: %.3f\n', ...
        regionName, iDir, peakLag, peakCorr);
end

fprintf('\n========== summary (cortex) ==========\n');
fprintf('peak lags (s):\n');  disp(peakLags);
fprintf('peak corrs   :\n');  disp(peakCorrs);
fprintf('lag cis (s)  :\n');  disp(lagCIAll);

savePath = fullfile('X:\David\AnalysesData\InterneuronAnalyses\Lab Meeting Pres\4 aninals run cross correlation, no chunking, pop-wise', ...
    'runCrossCorrelation_savedOutputs_all4Animals.mat');

save(savePath, ...
    'xcorrResults', ...
    'peakLags', ...
    'peakCorrs', ...
    'lagCIAll', ...
    'permLagCell', ...
    'animalLabels');

fprintf('\nsaved plotting-ready outputs to:\n%s\n', savePath);

end

% ================= helper =================
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
