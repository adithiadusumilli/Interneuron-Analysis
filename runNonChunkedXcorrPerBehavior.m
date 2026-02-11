function runNonChunkedXcorrPerBehavior(baseDirs, labelType, behaviors, maxLagSecs, numShifts, numPerms)
% non-chunked (full time series) cross-correlation between cortex interneuron and pyramidal population activity, computed separately within each behavior label type

% this is a hybrid of runCrossCorrelation (population mean int vs pyr from cortexFRs + AA_classifications) and chunked per-behavior tiled logic (one tile per behavior, save per behavior outputs)

% important: this is NOT event-centered like chunked, it j restricts to timepoints belonging to a behavior, concatenates those timepoints into vectors, then computes xc(lag) over those vectors

% label types:
%   - "umap": 1..7 (computed from regionAssignmentsFiltered mapped via origDownsampEMGInd)
%   - "manual": canonical 0..10 (computed from behvLabelsNoArt + analyzedBehaviors -> manBehvNames)
%   - "classifier": canonical 0..10 (computed from classifierLabels + classifierBehvs -> manBehvNames)

% inputs:
%   baseDirs : cell array of ProcessedData folders (one per animal/session)
%   labelType : "umap" | "manual" | "classifier"
%   behaviors : vector of behaviors to compute (default depends on labelType)
%   maxLagSecs : max lag in seconds (default 0.5)
%   numShifts : circular-shift null draws for correlation bounds (default 100)
%   numPerms : label-permutation draws for lag ci (default 100)

% outputs:
%   - per session: a 1 x nBehaviors tiled figure with xc(lag) + null bounds
%   - per session: saves a .mat containing per-behavior xc, peak lag, and null summaries

% example:
%   baseDirs = {
%     'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData', ...
%     'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData', ...
%     'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData', ...
%     'X:\David\ArenaRecordings\D043-020425-ArenaRecording\ProcessedData'
%   };
%   runNonChunkedXcorrPerBehavior(baseDirs, "umap", 1:7, 0.5, 100, 100);
%   runNonChunkedXcorrPerBehavior(baseDirs, "manual", 0:10, 0.5, 100, 100);
%   runNonChunkedXcorrPerBehavior(baseDirs, "classifier", 0:10, 0.5, 100, 100);

    if nargin < 2 || isempty(labelType),  labelType = "umap"; end
    labelType = lower(string(labelType));

    if nargin < 4 || isempty(maxLagSecs), maxLagSecs = 0.5; end
    if nargin < 5 || isempty(numShifts),  numShifts  = 100; end
    if nargin < 6 || isempty(numPerms),   numPerms   = 100; end

    if nargin < 3 || isempty(behaviors)
        if labelType == "umap"
            behaviors = 1:7;
        else
            behaviors = 0:10;
        end
    end

    binSize = 0.001;  % 1 ms bins
    nSess   = numel(baseDirs);
    nB      = numel(behaviors);

    % load centralized classifications file
    conslidatedDataFoler = 'X:\David\AnalysesData';
    load(fullfile(conslidatedDataFoler, 'AA_classifications.mat'), 'classifications');

    % known base folders (match AA_classifications row order)
    animalFolders = {
        'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData', ...
        'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData', ...
        'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData', ...
        'X:\David\ArenaRecordings\D043-020425-ArenaRecording\ProcessedData'
    };

    for iDir = 1:nSess
        baseDir = baseDirs{iDir};
        fprintf('\n=== processing session %d/%d: %s ===\n', iDir, nSess, baseDir);

        frFile = fullfile(baseDir, 'NeuralFiringRates1msBins10msGauss.mat');
        if ~isfile(frFile)
            warning('missing firing rate file in %s. skipping.', baseDir);
            continue;
        end

        % load cortex firing rates + inds
        F = load(frFile, 'cortexFRs', 'cortexInds');
        cortexFRs  = F.cortexFRs;
        cortexInds = F.cortexInds(:);

        frMatrix = cortexFRs; % neurons x time
        regionInds = cortexInds; % neuronDataStruct indices for cortex units
        Tfr = size(frMatrix, 2);

        % match classifications row
        matchRow = find(contains(animalFolders, baseDir), 1);
        if isempty(matchRow)
            warning('could not match baseDir to animalFolders list. skipping.');
            continue;
        end

        % cortex is iRegion = 1
        iRegion = 1;
        neuronType = classifications{matchRow, iRegion}; % 1=int, 0=pyr
        if isempty(neuronType)
            warning('no classification data for %s. skipping.', baseDir);
            continue;
        end

        regionClass = neuronType(regionInds);

        intFRs = frMatrix(regionClass == 1, :);
        pyrFRs = frMatrix(regionClass == 0, :);

        if isempty(intFRs) || isempty(pyrFRs)
            warning('no valid int/pyr cortex units for %s. skipping.', baseDir);
            continue;
        end

        % population means across neurons (this is the "population analysis")
        meanIntRaw = mean(intFRs, 1, 'omitnan');
        meanPyrRaw = mean(pyrFRs, 1, 'omitnan');

        % build full-length behavior vector aligned to fr time base
        % note: this uses origDownsampEMGInd mapping logic like your canonical mapper,
        % but now for timepoints (not transitions)
        [labels1k, meta] = buildTimepointBehaviorLabels(baseDir, Tfr, labelType);

        % make tiled figure
        figure('Color','w', 'Name', sprintf('non-chunked xcorr by %s | sess %d', labelType, iDir));
        tiledlayout(1, nB, 'TileSpacing','compact','Padding','compact');

        % outputs per session
        out = struct();
        out.baseDir = baseDir;
        out.labelType = labelType;
        out.behaviors = behaviors;
        out.binSize = binSize;
        out.maxLagSecs = maxLagSecs;
        out.numShifts = numShifts;
        out.numPerms = numPerms;
        out.meta = meta;

        out.beh = repmat(struct('beh',[], 'nTimepoints',[], 'timeIdx',[], 'lagsSec',[], 'xc',[], 'peakLagSec',[], 'peakCorr',[], 'ctrlCorrs',[], 'ctrlCorrCI',[], 'permPeakLags',[], 'lagCI',[]), 1, nB);

        for bIdx = 1:nB
            beh = behaviors(bIdx);

            nexttile; hold on;

            timeIdx = find(labels1k == beh);
            if isempty(timeIdx)
                title(sprintf('beh %d (n=0)', beh));
                axis off;
                out.beh(bIdx).beh = beh;
                out.beh(bIdx).nTimepoints = 0;
                continue;
            end

            % concatenate behavior timepoints into vectors
            intVec = meanIntRaw(timeIdx);
            pyrVec = meanPyrRaw(timeIdx);

            % compute xc(lag)
            [lagsSec, xc, peakLag, peakCorr] = computeManualXCorrVec(intVec, pyrVec, binSize, maxLagSecs);

            % control distribution for correlation bounds (circular shifts)
            ctrlCorrs = nan(1, numShifts);
            minShiftBins = round(30 / binSize); % 30 sec minimum shift
            maxShiftBinsOk = numel(intVec) - minShiftBins;

            if maxShiftBinsOk >= minShiftBins
                for s = 1:numShifts
                    shiftAmt = randi([minShiftBins, maxShiftBinsOk]);
                    intShift = circshift(intVec, shiftAmt);
                    v = ~isnan(intShift) & ~isnan(pyrVec);
                    if nnz(v) > 10
                        ctrlCorrs(s) = corr(intShift(v)', pyrVec(v)');
                    end
                end
                ctrlCI = prctile(ctrlCorrs, [2.5 97.5]);
            else
                ctrlCI = [NaN NaN];
            end

            % permutation distribution for lag (shuffle neuron labels within cortex)
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

                [~, ~, permPeakLag, ~] = computeManualXCorrVec(permIntVec, permPyrVec, binSize, maxLagSecs);
                permPeakLags(p) = permPeakLag;
            end

            goodPerm = ~isnan(permPeakLags);
            if any(goodPerm)
                lagCI = prctile(permPeakLags(goodPerm), [2.5 97.5]);
            else
                lagCI = [NaN NaN];
            end

            % -------- plot tile --------
            plot(lagsSec, xc, 'k', 'LineWidth', 2);
            xline(0, 'k:');

            if ~any(isnan(ctrlCI))
                yline(ctrlCI(1), '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);
                yline(ctrlCI(2), '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);
            end

            if ~isnan(peakLag)
                xline(peakLag, 'r--', 'LineWidth', 1.5);
            end

            if ~any(isnan(lagCI))
                xline(lagCI(1), '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);
                xline(lagCI(2), '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);
            end

            title(sprintf('beh %d (n=%d)', beh, numel(timeIdx)));
            xlabel('lag (s)');
            if bIdx == 1, ylabel('correlation'); end
            box off;

            % store outputs
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

            % save a per-behavior mat like chunked code
            outFileBeh = fullfile(baseDir, sprintf('nonchunked_xcorr_withNull_%s_beh%02d.mat', labelType, beh));
            lagsBins = round(lagsSec / binSize);
            nTimepoints = numel(timeIdx);
            save(outFileBeh, 'beh','labelType','nTimepoints','timeIdx','lagsBins','lagsSec','binSize', 'xc','peakLag','peakCorr','ctrlCorrs','ctrlCI','permPeakLags','lagCI');
        end

        sgtitle(sprintf('non-chunked xcorr by %s | maxLag=\\pm%.3f s | %s', labelType, maxLagSecs, shortTag(baseDir)));

        % save full session output
        outFile = fullfile(baseDir, sprintf('nonchunked_xcorr_by_%s_cortex_maxlag%.0fms.mat', labelType, 1000*maxLagSecs));
        save(outFile, 'out');
        fprintf('saved session output: %s\n', outFile);
    end
end

% =====================================================================
% helper: build full-length timepoint behavior labels aligned to fr bins
% =====================================================================
function [labels1k, meta] = buildTimepointBehaviorLabels(baseDir, Tfr, labelType)
% outputs:
%   labels1k: 1 x Tfr vector with labels (umap:1..7, manual/classifier:0..10), nan if unknown
%   meta: struct with mapping info and canonical names (for manual/classifier)

    labelType = lower(string(labelType));
    meta = struct();
    meta.labelType = labelType;

    labels1k = nan(1, Tfr);

    umapFile = fullfile(baseDir, 'UMAP.mat');
    if ~isfile(umapFile)
        warning('UMAP.mat not found in %s', baseDir);
        return;
    end

    U = load(umapFile, ...
        'regionAssignmentsFiltered', 'regionAssignmentFiltered', ...
        'behvLabelsNoArt', 'origDownsampEMGInd', ...
        'classifierLabels', 'classifierBehvs', ...
        'analyzedBehaviors');

    if ~isfield(U, 'origDownsampEMGInd') || isempty(U.origDownsampEMGInd)
        warning('origDownsampEMGInd missing/empty in UMAP.mat; cannot align labels.');
        return;
    end

    origInd = U.origDownsampEMGInd(:);

    % canonical behavior name list (shared across animals)
    manBehvNames = {'climbdown','climbup','eating','grooming', ...
                    'jumpdown','jumping','rearing','still','walkflat','walkgrid'};
    meta.manBehvNames = manBehvNames;

    % ----- choose raw label vector (reduced time base) -----
    labelVecReduced = [];

    switch labelType
        case "umap"
            if isfield(U,'regionAssignmentsFiltered') && ~isempty(U.regionAssignmentsFiltered)
                regionAssignmentsFiltered = U.regionAssignmentsFiltered;
            elseif isfield(U,'regionAssignmentFiltered') && ~isempty(U.regionAssignmentFiltered)
                regionAssignmentsFiltered = U.regionAssignmentFiltered;
            else
                warning('regionAssignment(s)Filtered missing; cannot build umap labels.');
                return;
            end

            regionAssignmentsFiltered = regionAssignmentsFiltered(:);

            % map raw region codes -> contiguous 1..7 like your mapper
            regionCodes = unique(regionAssignmentsFiltered(~isnan(regionAssignmentsFiltered)));
            regionCodes = regionCodes(:);
            if numel(regionCodes) ~= 7
                warning('expected 7 umap regions, found %d. still proceeding with contiguous mapping.', numel(regionCodes));
            end
            meta.regionCodes = regionCodes;

            labelVecReduced = nan(size(regionAssignmentsFiltered));
            for i = 1:numel(regionAssignmentsFiltered)
                rv = regionAssignmentsFiltered(i);
                if isnan(rv), continue; end
                idx = find(regionCodes == rv, 1);
                if ~isempty(idx)
                    labelVecReduced(i) = idx;
                end
            end

        case "manual"
            if ~isfield(U,'behvLabelsNoArt') || isempty(U.behvLabelsNoArt)
                warning('behvLabelsNoArt missing; cannot build manual labels.');
                return;
            end
            if ~isfield(U,'analyzedBehaviors') || isempty(U.analyzedBehaviors)
                warning('analyzedBehaviors missing; cannot remap manual labels to canonical.');
                return;
            end

            behvLabelsNoArt = U.behvLabelsNoArt(:);
            analyzedBehaviors = U.analyzedBehaviors;

            % build remap: per-animal manual index -> canonical index (0..10)
            nManualBehv = numel(analyzedBehaviors);
            manBehvNumbers = zeros(1, nManualBehv); % 0 = not tracked in canonical list

            for iBehv = 1:nManualBehv
                thisName = analyzedBehaviors{iBehv};
                idx = find(strcmp(thisName, manBehvNames), 1);
                if isempty(idx)
                    manBehvNumbers(iBehv) = 0;
                else
                    manBehvNumbers(iBehv) = idx;
                end
            end
            meta.analyzedBehaviors = analyzedBehaviors;
            meta.manualRemap = manBehvNumbers;

            % apply remap to reduced label vector
            labelVecReduced = nan(size(behvLabelsNoArt));
            for i = 1:numel(behvLabelsNoArt)
                v = behvLabelsNoArt(i);
                if isnan(v)
                    continue;
                elseif v == 0
                    labelVecReduced(i) = 0;
                else
                    if v >= 1 && v <= numel(manBehvNumbers)
                        labelVecReduced(i) = manBehvNumbers(v);
                    else
                        labelVecReduced(i) = nan;
                    end
                end
            end

        case "classifier"
            if ~isfield(U,'classifierLabels') || isempty(U.classifierLabels)
                warning('classifierLabels missing; cannot build classifier labels.');
                return;
            end
            if ~isfield(U,'classifierBehvs') || isempty(U.classifierBehvs)
                warning('classifierBehvs missing; cannot remap classifier labels to canonical.');
                return;
            end

            classifierLabels = U.classifierLabels(:);
            classifierBehvs = U.classifierBehvs;

            % build remap: per-animal classifier index -> canonical index (0..10)
            nClassBehv = numel(classifierBehvs);
            classBehvNumbers = zeros(1, nClassBehv);

            for iBehv = 1:nClassBehv
                thisName = classifierBehvs{iBehv};
                idx = find(strcmp(thisName, manBehvNames), 1);
                if isempty(idx)
                    classBehvNumbers(iBehv) = 0;
                else
                    classBehvNumbers(iBehv) = idx;
                end
            end
            meta.classifierBehvs = classifierBehvs;
            meta.classifierRemap = classBehvNumbers;

            % apply remap to reduced label vector
            labelVecReduced = nan(size(classifierLabels));
            for i = 1:numel(classifierLabels)
                v = classifierLabels(i);
                if isnan(v)
                    continue;
                elseif v == 0
                    labelVecReduced(i) = 0;
                else
                    if v >= 1 && v <= numel(classBehvNumbers)
                        labelVecReduced(i) = classBehvNumbers(v);
                    else
                        labelVecReduced(i) = nan;
                    end
                end
            end

        otherwise
            error('labelType must be "umap", "manual", or "classifier".');
    end

    % ----- align reduced labels into full fr timeline using origDownsampEMGInd -----
    n = min(numel(origInd), numel(labelVecReduced));
    origInd = origInd(1:n);
    labelVecReduced = labelVecReduced(1:n);

    ok = origInd >= 1 & origInd <= Tfr & ~isnan(labelVecReduced);
    labels1k(origInd(ok)) = labelVecReduced(ok);

    meta.nReduced = numel(labelVecReduced);
    meta.nAligned = nnz(ok);
end

% =====================================================================
% helper: compute xc for two vectors by scanning lags (removes nans)
% =====================================================================
function [lagsSec, xc, peakLagSec, peakCorr] = computeManualXCorrVec(intVec, pyrVec, binSize, maxLagSecs)
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
            xc(li) = corr(intSeg(v), pyrSeg(v));
        end
    end

    [peakCorr, peakIdx] = max(xc);
    if isempty(peakIdx) || isnan(peakCorr)
        peakLagSec = NaN;
    else
        peakLagSec = lagsSec(peakIdx);
    end
end

% =====================================================================
% helper: short label for sgtitle
% =====================================================================
function tag = shortTag(baseDir)
    [~, folderName] = fileparts(baseDir);
    m = regexp(folderName, 'D\d+', 'match', 'once');
    if isempty(m), tag = folderName; else, tag = m; end
end
