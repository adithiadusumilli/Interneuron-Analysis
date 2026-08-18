function runNonChunkedXcorrPerBehavior_saveOnly_getMouseDataNames(mouseIDs, baseSessionNames, probeRegions, labelType, behaviors, maxLagSecs, numShifts, numPerms)
% non-chunked (full time series) population-averaged cross-correlation
% computed separately per behv using David's getMouseDataNames

% IMPORTANT:
%   no trial averaging because the analysis is not event-centered.
%   for each behavior, the code:
%       1. averages across interneurons -> 1 full-session int population trace
%       2. averages across pyramidal neurons -> 1 full-session pyr population trace
%       3. selects the timepoints labeled as that behavior
%       4. computes xcorr on those behavior-specific vectors


if nargin < 4 || isempty(labelType), labelType = "classifier"; end
labelType = lower(string(labelType));
if labelType ~= "classifier"
    error('This updated version only supports labelType = "classifier".');
end

if nargin < 5 || isempty(behaviors), behaviors = 0:10; end
if nargin < 6 || isempty(maxLagSecs), maxLagSecs = 0.5; end
if nargin < 7 || isempty(numShifts), numShifts = 100; end
if nargin < 8 || isempty(numPerms), numPerms = 100; end

nSess = numel(mouseIDs);
if numel(baseSessionNames) ~= nSess || numel(probeRegions) ~= nSess
    error('mouseIDs, baseSessionNames, and probeRegions must have the same length.');
end

binSize = 0.001;
nB = numel(behaviors);

consolidatedDataFolder = 'X:\David\AnalysesData';
classificationFile = fullfile(consolidatedDataFolder, 'AA_classifications.mat');
if ~isfile(classificationFile)
    error('Missing classifications file: %s', classificationFile);
end
load(classificationFile, 'classifications');

results = struct();
results.mouseIDs = mouseIDs;
results.baseSessionNames = baseSessionNames;
results.probeRegions = probeRegions;
results.labelType = labelType;
results.behaviors = behaviors;
results.binSize = binSize;
results.maxLagSecs = maxLagSecs;
results.numShifts = numShifts;
results.numPerms = numPerms;
results.sessions = repmat(struct( ...
    'mouseID', '', ...
    'baseSessionName', '', ...
    'probeRegion', '', ...
    'baseDir', '', ...
    'processedDataFolder', '', ...
    'sessionTag', '', ...
    'meta', struct(), ...
    'beh', []), 1, nSess);

for iDir = 1:nSess

    mouseID = mouseIDs{iDir};
    baseSessionName = baseSessionNames{iDir};
    probeRegion = probeRegions{iDir};

    fprintf('\n============================================================\n');
    fprintf('processing session %d/%d: %s\n', iDir, nSess, mouseID);
    fprintf('baseSessionName: %s\n', baseSessionName);
    fprintf('probeRegion: %s\n', probeRegion);
    fprintf('============================================================\n');

    dataNames = getMouseDataNames(mouseID, baseSessionName, probeRegion);

    if ~isfield(dataNames, 'processedDataFolder') || isempty(dataNames.processedDataFolder)
        warning('getMouseDataNames returned no processedDataFolder for %s. Skipping.', mouseID);
        continue;
    end

    baseDir = dataNames.processedDataFolder;

    if isfield(dataNames, 'NeuralFiringRates1msBins10msGauss') && ...
            ~isempty(dataNames.NeuralFiringRates1msBins10msGauss)
        frFile = dataNames.NeuralFiringRates1msBins10msGauss;
    else
        frFile = fullfile(baseDir, 'NeuralFiringRates1msBins10msGauss.mat');
    end

    fprintf('processedDataFolder:\n%s\n', baseDir);

    if ~isfile(frFile)
        warning('missing firing rate file for %s: %s. Skipping.', mouseID, frFile);
        continue;
    end

    F = load(frFile, 'cortexFRs', 'cortexInds');
    cortexFRs = F.cortexFRs;
    cortexInds = F.cortexInds(:);

    frMatrix = cortexFRs;
    regionInds = cortexInds;
    Tfr = size(frMatrix, 2);

    matchRow = getClassificationRow(mouseID);

    if matchRow > size(classifications, 1)
        warning('classification row %d for %s is unavailable. Skipping.', matchRow, mouseID);
        continue;
    end

    neuronType = classifications{matchRow, 1};
    if isempty(neuronType)
        warning('no classification data for %s. Skipping.', mouseID);
        continue;
    end

    if max(regionInds) > numel(neuronType)
        warning('cortexInds exceed classification vector length for %s. Skipping.', mouseID);
        continue;
    end

    regionClass = neuronType(regionInds);

    intFRs = frMatrix(regionClass == 1, :);
    pyrFRs = frMatrix(regionClass == 0, :);

    if isempty(intFRs) || isempty(pyrFRs)
        warning('no valid int/pyr cortex units for %s. Skipping.', mouseID);
        continue;
    end

    fprintf('nInt=%d | nPyr=%d | total cortex=%d\n', ...
        size(intFRs,1), size(pyrFRs,1), numel(regionClass));

    meanIntRaw = mean(intFRs, 1, 'omitnan');
    meanPyrRaw = mean(pyrFRs, 1, 'omitnan');

    [labels1k, meta] = buildTimepointBehaviorLabels(baseDir, Tfr, labelType);

    out = struct();
    out.mouseID = mouseID;
    out.baseSessionName = baseSessionName;
    out.probeRegion = probeRegion;
    out.baseDir = baseDir;
    out.labelType = labelType;
    out.behaviors = behaviors;
    out.binSize = binSize;
    out.maxLagSecs = maxLagSecs;
    out.numShifts = numShifts;
    out.numPerms = numPerms;
    out.meta = meta;

    out.beh = repmat(struct( ...
        'beh', [], ...
        'nTimepoints', [], ...
        'timeIdx', [], ...
        'lagsSec', [], ...
        'xc', [], ...
        'peakLagSec', [], ...
        'peakCorr', [], ...
        'ctrlCorrs', [], ...
        'ctrlCorrCI', [], ...
        'permPeakLags', [], ...
        'lagCI', []), 1, nB);

    for bIdx = 1:nB

        beh = behaviors(bIdx);
        timeIdx = find(labels1k == beh);

        if isempty(timeIdx)
            out.beh(bIdx).beh = beh;
            out.beh(bIdx).nTimepoints = 0;
            continue;
        end

        intVec = meanIntRaw(timeIdx);
        pyrVec = meanPyrRaw(timeIdx);

        [lagsSec, xc, peakLag, peakCorr] = computeManualXCorrVec( ...
            intVec, pyrVec, binSize, maxLagSecs);

        ctrlCorrs = nan(1, numShifts);
        minShiftBins = round(30 / binSize);
        maxShiftBinsOk = numel(meanIntRaw) - minShiftBins;

        if maxShiftBinsOk >= minShiftBins
            for s = 1:numShifts
                shiftAmt = randi([minShiftBins, maxShiftBinsOk]);

                intShiftFull = circshift(meanIntRaw, shiftAmt);
                intShiftVec = intShiftFull(timeIdx);

                v = ~isnan(intShiftVec) & ~isnan(pyrVec);
                if nnz(v) > 10
                    x = intShiftVec(v);
                    y = pyrVec(v);
                    ctrlCorrs(s) = corr(x(:), y(:));
                end
            end

            goodCtrl = isfinite(ctrlCorrs);
            if any(goodCtrl)
                ctrlCI = prctile(ctrlCorrs(goodCtrl), [2.5 97.5]);
            else
                ctrlCI = [NaN NaN];
            end
        else
            ctrlCI = [NaN NaN];
        end

        permPeakLags = nan(1, numPerms);

        for p = 1:numPerms
            permLabels = regionClass(randperm(numel(regionClass)));

            permIntFRs = frMatrix(permLabels == 1, :);
            permPyrFRs = frMatrix(permLabels == 0, :);

            if isempty(permIntFRs) || isempty(permPyrFRs)
                continue;
            end

            permMeanInt = mean(permIntFRs, 1, 'omitnan');
            permMeanPyr = mean(permPyrFRs, 1, 'omitnan');

            permIntVec = permMeanInt(timeIdx);
            permPyrVec = permMeanPyr(timeIdx);

            [~, ~, permPeakLag, ~] = computeManualXCorrVec( ...
                permIntVec, permPyrVec, binSize, maxLagSecs);

            permPeakLags(p) = permPeakLag;
        end

        goodPerm = isfinite(permPeakLags);
        if any(goodPerm)
            lagCI = prctile(permPeakLags(goodPerm), [2.5 97.5]);
        else
            lagCI = [NaN NaN];
        end

        out.beh(bIdx).beh = beh;
        out.beh(bIdx).nTimepoints = numel(timeIdx);
        out.beh(bIdx).timeIdx = timeIdx;
        out.beh(bIdx).lagsSec = lagsSec;
        out.beh(bIdx).xc = xc;
        out.beh(bIdx).peakLagSec = peakLag;
        out.beh(bIdx).peakCorr = peakCorr;
        out.beh(bIdx).ctrlCorrs = ctrlCorrs;
        out.beh(bIdx).ctrlCorrCI = ctrlCI;
        out.beh(bIdx).permPeakLags = permPeakLags(goodPerm);
        out.beh(bIdx).lagCI = lagCI;

        fprintf('  behavior %2d | nTimepoints=%d | peakLag=%+.4f s | peakCorr=%.4f | validPerm=%d/%d\n', ...
            beh, numel(timeIdx), peakLag, peakCorr, nnz(goodPerm), numPerms);
    end

    results.sessions(iDir).mouseID = mouseID;
    results.sessions(iDir).baseSessionName = baseSessionName;
    results.sessions(iDir).probeRegion = probeRegion;
    results.sessions(iDir).baseDir = baseDir;
    results.sessions(iDir).processedDataFolder = baseDir;
    results.sessions(iDir).sessionTag = mouseID;
    results.sessions(iDir).meta = meta;
    results.sessions(iDir).beh = out.beh;
end

savePath = fullfile( ...
    consolidatedDataFolder, ...
    sprintf('nonchunked_xcorr_by_%s_cortex_all6Sessions_saved.mat', labelType));

save(savePath, 'results', '-v7.3');
fprintf('\nsaved plotting-ready outputs:\n%s\n', savePath);

end

function [labels1k, meta] = buildTimepointBehaviorLabels(baseDir, Tfr, labelType)

labelType = lower(string(labelType));
meta = struct();
meta.labelType = labelType;
labels1k = nan(1, Tfr);

if labelType ~= "classifier"
    error('buildTimepointBehaviorLabels only supports "classifier".');
end

umapFile = fullfile(baseDir, 'UMAP.mat');

if ~isfile(umapFile)
    warning('UMAP.mat not found in %s', baseDir);
    return;
end

U = load(umapFile, ...
    'origDownsampEMGInd', ...
    'classifierLabels', ...
    'classifierBehvs');

if ~isfield(U, 'origDownsampEMGInd') || isempty(U.origDownsampEMGInd)
    warning('origDownsampEMGInd missing/empty in UMAP.mat; cannot align labels.');
    return;
end

if ~isfield(U, 'classifierLabels') || isempty(U.classifierLabels)
    warning('classifierLabels missing; cannot build classifier labels.');
    return;
end

if ~isfield(U, 'classifierBehvs') || isempty(U.classifierBehvs)
    warning('classifierBehvs missing; cannot remap classifier labels to canonical.');
    return;
end

origInd = U.origDownsampEMGInd(:);
classifierLabels = U.classifierLabels(:);
classifierBehvs = U.classifierBehvs;

manBehvNames = {'climbdown','climbup','eating','grooming', ...
                'jumpdown','jumping','rearing','still','walkflat','walkgrid'};

meta.manBehvNames = manBehvNames;

canonicalLookup = containers.Map;
canonicalLookup('climbdown') = 1;
canonicalLookup('climbup') = 2;
canonicalLookup('eating') = 3;
canonicalLookup('eat') = 3;
canonicalLookup('grooming') = 4;
canonicalLookup('groom') = 4;
canonicalLookup('jumpdown') = 5;
canonicalLookup('jumping') = 6;
canonicalLookup('jumpacross') = 6;
canonicalLookup('rearing') = 7;
canonicalLookup('rear') = 7;
canonicalLookup('still') = 8;
canonicalLookup('walkflat') = 9;
canonicalLookup('walkgrid') = 10;

nClassBehv = numel(classifierBehvs);
classBehvNumbers = zeros(1, nClassBehv);

for iBehv = 1:nClassBehv
    thisName = classifierBehvs{iBehv};
    cleanName = lower(strrep(strrep(thisName, ' ', ''), '_', ''));

    if isKey(canonicalLookup, cleanName)
        classBehvNumbers(iBehv) = canonicalLookup(cleanName);
    else
        classBehvNumbers(iBehv) = 0;
    end
end

meta.classifierBehvs = classifierBehvs;
meta.classifierRemap = classBehvNumbers;

labelVecReduced = nan(size(classifierLabels));

for i = 1:numel(classifierLabels)
    v = classifierLabels(i);

    if isnan(v)
        continue;
    elseif v == 0
        labelVecReduced(i) = 0;
    elseif v >= 1 && v <= numel(classBehvNumbers)
        labelVecReduced(i) = classBehvNumbers(v);
    else
        labelVecReduced(i) = NaN;
    end
end

n = min(numel(origInd), numel(labelVecReduced));
origInd = origInd(1:n);
labelVecReduced = labelVecReduced(1:n);

ok = origInd >= 1 & origInd <= Tfr & ~isnan(labelVecReduced);
labels1k(origInd(ok)) = labelVecReduced(ok);

meta.nReduced = numel(labelVecReduced);
meta.nAligned = nnz(ok);

end

function [lagsSec, xc, peakLagSec, peakCorr] = computeManualXCorrVec( ...
    intVec, pyrVec, binSize, maxLagSecs)

maxLagBins = round(maxLagSecs / binSize);
lagsBins = -maxLagBins:maxLagBins;
lagsSec = lagsBins * binSize;
xc = nan(size(lagsBins));

intVec = intVec(:)';
pyrVec = pyrVec(:)';

for li = 1:numel(lagsBins)
    lag = lagsBins(li);

    if lag < 0
        intSeg = intVec(1:end+lag);
        pyrSeg = pyrVec(1-lag:end);
    elseif lag > 0
        intSeg = intVec(1+lag:end);
        pyrSeg = pyrVec(1:end-lag);
    else
        intSeg = intVec;
        pyrSeg = pyrVec;
    end

    v = ~isnan(intSeg) & ~isnan(pyrSeg);

    if nnz(v) > 10
        x = intSeg(v);
        y = pyrSeg(v);
        xc(li) = corr(x(:), y(:));
    end
end

[peakCorr, peakIdx] = max(xc);

if isempty(peakIdx) || isnan(peakCorr)
    peakLagSec = NaN;
else
    peakLagSec = lagsSec(peakIdx);
end

end

function matchRow = getClassificationRow(mouseID)

switch upper(string(mouseID))
    case "D026"
        matchRow = 1;
    case "D020"
        matchRow = 2;
    case "D024"
        matchRow = 3;
    case "D043"
        matchRow = 4;
    case "D050"
        matchRow = 5;
    case "D054"
        matchRow = 6;
    otherwise
        error('Unknown mouseID for AA_classifications ordering: %s', mouseID);
end

end
