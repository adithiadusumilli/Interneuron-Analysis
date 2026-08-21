% chunked TRIAL-AVERAGED population cross-correlation per behavior
% UPDATED TO MATCH REGULAR CHUNKED POPAVG PIPELINE
%
% jobInd:
%   0       real
%   1:300   random-label permutations
%   301:400 shifted controls (shift 1:100)
%
% The 100 shifted controls are pooled separately for EACH animal x behavior
% to define one 95th-percentile zero-lag correlation threshold. Every
% permutation for that same animal x behavior is redrawn until its peak
% correlation reaches that threshold or maxPermTries is reached.
%
% Bayes-ready permutation blocks:
%   1:100   H0
%   101:200 H+50
%   201:300 H-50
%
% IMPORTANT: run shift jobs 301:400 BEFORE permutation jobs 1:300.

clc;

%% ---------------- settings ----------------
baseDirs = {
    '/home/asa7288/Transfer/D026', ...
    '/home/asa7288/Transfer/D020', ...
    '/home/asa7288/Transfer/D024', ...
    '/home/asa7288/Transfer/D043', ...
    '/home/asa7288/Transfer/D050', ...
    '/home/asa7288/Transfer/D054'
};

labelsToUse = 'classifier';
labelTag = lower(labelsToUse);

rng(jobInd)

chunkHalf = 200;
maxLagSecs = 0.5;
doBaselineNorm = true;
maxPermTries = 1000;
maxShiftRetries = 5000;

switch lower(labelsToUse)
    case 'umap'
        behaviors = 1:7;
    case 'manual'
        behaviors = 1:10;
    case 'classifier'
        behaviors = 1:10;
    otherwise
        error('unknown labelsToUse: %s', labelsToUse);
end

%% ---------------- job type ----------------
if jobInd == 0
    jobType = "real";
    permInd = 0;
    shiftInd = 0;
elseif jobInd >= 1 && jobInd <= 300
    jobType = "perm";
    permInd = jobInd;
    shiftInd = 0;
elseif jobInd >= 301 && jobInd <= 400
    jobType = "shift";
    permInd = 0;
    shiftInd = jobInd - 300;
else
    error('jobInd must be 0, 1:300, or 301:400');
end

%% ---------------- sessions ----------------
for iDir = 1:numel(baseDirs)
    baseDir = baseDirs{iDir};
    outDir = fullfile(baseDir,'quest_runs');
    if ~exist(outDir,'dir'), mkdir(outDir); end

    fprintf('\n=== %s | job %d | %s ===\n',baseDir,jobInd,jobType);

    if jobType == "shift"
        S = load(fullfile(baseDir,'EMG_Neural_AllChannels.mat'), ...
            'pyrCxWinCell','intCxWinCell','tAxis', ...
            'validTransitionsCell','validTransitionsNeurShiftedCell');

        syncData = load(fullfile(baseDir,'VideoSyncFrames.mat'), ...
            'frameEMGSamples','frameNeuropixelSamples');
        frameEMGSamples = syncData.frameEMGSamples;
        frameNeuropixelSamples = syncData.frameNeuropixelSamples;

        emgData = load(fullfile(baseDir,'EMG1ms.mat'),'downsampEMG');
        totalLength = size(emgData.downsampEMG,2);
        minShift = 30000;
    else
        S = load(fullfile(baseDir,'EMG_Neural_AllChannels.mat'), ...
            'pyrCxWinCell','intCxWinCell','tAxis');
    end

    L = load(fullfile(baseDir,'transitionCanonicalBehaviorLabels.mat'), ...
        'regionLabelsPerTransition','manualLabelsPerTransition','classifierLabelsPerTransition');

    switch lower(labelsToUse)
        case 'umap'
            labelsPerTransition = L.regionLabelsPerTransition;
        case 'manual'
            labelsPerTransition = L.manualLabelsPerTransition;
        case 'classifier'
            labelsPerTransition = L.classifierLabelsPerTransition;
    end

    % Upstream EMG_Neural_AllChannels.mat contains only good channels.
    channelsToUseThisSess = 1:numel(S.pyrCxWinCell);

    tAxis = S.tAxis(:)';
    binSize = 0.001;
    T = numel(tAxis);
    [~,cIdx] = min(abs(tAxis));
    pyrStart = cIdx - chunkHalf;
    pyrEnd = cIdx + chunkHalf;

    reqMaxLagBins = round(maxLagSecs/binSize);
    maxLagLeft = pyrStart - 1;
    maxLagRight = T - pyrEnd;
    maxLagBins = min([reqMaxLagBins,maxLagLeft,maxLagRight]);
    lags = -maxLagBins:maxLagBins;

    for bIdx = 1:numel(behaviors)
        beh = behaviors(bIdx);
        fprintf('  behavior %d\n',beh);

        %% ====================== SHIFT ======================
        if jobType == "shift"
            [cortexIntFR,okI] = loadCortexIntFR(baseDir);
            if ~okI
                warning('failed to load cortex interneuron FR for %s; beh %d skipped',baseDir,beh);
                continue;
            end

            requestedShiftCol = shiftInd + 1;
            usedShiftCol = NaN;
            usedShiftAmtMs = NaN;
            usedFallbackShift = false;
            foundValidShift = false;
            retryCount = 0;

            [pyr_byCh_try,int_byCh_try,nTrials_try] = buildShiftedBehaviorSample( ...
                S,labelsPerTransition,beh,requestedShiftCol,cortexIntFR,tAxis, ...
                doBaselineNorm,false,[],[],[],[],[],channelsToUseThisSess);

            if nTrials_try > 0
                [pyrAvgTrace_try,intAvgTrace_try,nTrialsUsed_try] = ...
                    averageAcrossTrialsAndChannels(pyr_byCh_try,int_byCh_try);
                if nTrialsUsed_try > 0 && ~isempty(pyrAvgTrace_try) && ~isempty(intAvgTrace_try)
                    xcZeroLag_try = zeroLagTrialAverage(pyrAvgTrace_try,intAvgTrace_try,pyrStart,pyrEnd);
                else
                    xcZeroLag_try = NaN;
                end
                if isfinite(xcZeroLag_try)
                    pyrAvgTrace = pyrAvgTrace_try;
                    intAvgTrace = intAvgTrace_try;
                    nTrialsThisBehavior = nTrials_try;
                    nTrialsUsed = nTrialsUsed_try;
                    xcZeroLag = xcZeroLag_try;
                    usedShiftCol = requestedShiftCol;
                    foundValidShift = true;
                end
            end

            if ~foundValidShift
                for retryIdx = 1:maxShiftRetries
                    retryCount = retryIdx;
                    randShift = randi([minShift,totalLength-minShift]);
                    if rand() > 0.5, randShift = -randShift; end

                    [pyr_byCh_try,int_byCh_try,nTrials_try] = buildShiftedBehaviorSample( ...
                        S,labelsPerTransition,beh,NaN,cortexIntFR,tAxis, ...
                        doBaselineNorm,true,randShift,frameEMGSamples, ...
                        frameNeuropixelSamples,totalLength,baseDir,channelsToUseThisSess);

                    if nTrials_try < 1, continue; end
                    [pyrAvgTrace_try,intAvgTrace_try,nTrialsUsed_try] = ...
                        averageAcrossTrialsAndChannels(pyr_byCh_try,int_byCh_try);
                    if nTrialsUsed_try > 0 && ~isempty(pyrAvgTrace_try) && ~isempty(intAvgTrace_try)
                        xcZeroLag_try = zeroLagTrialAverage(pyrAvgTrace_try,intAvgTrace_try,pyrStart,pyrEnd);
                    else
                        xcZeroLag_try = NaN;
                    end
                    if isfinite(xcZeroLag_try)
                        pyrAvgTrace = pyrAvgTrace_try;
                        intAvgTrace = intAvgTrace_try;
                        nTrialsThisBehavior = nTrials_try;
                        nTrialsUsed = nTrialsUsed_try;
                        xcZeroLag = xcZeroLag_try;
                        usedShiftAmtMs = randShift;
                        usedFallbackShift = true;
                        foundValidShift = true;
                        break;
                    end
                end
            end

            if ~foundValidShift
                warning('no valid shifted sample for %s beh %d',baseDir,beh);
                continue;
            end

            runtimeSec = NaN;
            outFile = fullfile(outDir,sprintf( ...
                'concatCrossCorrPerCanonicalBehavior_%s_shift_%03d_zerolag_beh%02d.mat', ...
                labelTag,shiftInd,beh));

            save(outFile,'binSize','xcZeroLag','chunkHalf','doBaselineNorm', ...
                'jobInd','permInd','shiftInd','jobType','baseDir','runtimeSec', ...
                'beh','nTrialsThisBehavior','nTrialsUsed','labelsToUse', ...
                'requestedShiftCol','usedShiftCol','usedShiftAmtMs', ...
                'usedFallbackShift','maxShiftRetries','retryCount', ...
                'pyrAvgTrace','intAvgTrace','channelsToUseThisSess');

            fprintf('    shift %d beh %d xcZeroLag=%.6f saved\n',shiftInd,beh,xcZeroLag);
            continue;
        end

        %% ====================== REAL ======================
        if jobType == "real"
            [pyr_byCh,int_byCh,nTrialsThisBehavior] = buildRealOrPermBehaviorSample( ...
                S,labelsPerTransition,beh,channelsToUseThisSess,tAxis,doBaselineNorm,false);

            if nTrialsThisBehavior < 1, continue; end
            [pyrAvgTrace,intAvgTrace,nTrialsUsed] = averageAcrossTrialsAndChannels(pyr_byCh,int_byCh);
            if nTrialsUsed < 1, continue; end

            tic;
            [xc,peakLagSec,peakCorr] = lagSweepTrialAverage( ...
                pyrAvgTrace,intAvgTrace,lags,pyrStart,pyrEnd,T,binSize);
            runtimeSec = toc;

            outFile = fullfile(outDir,sprintf( ...
                'concatCrossCorrPerCanonicalBehavior_%s_unperm_beh%02d.mat',labelTag,beh));
            save(outFile,'lags','binSize','xc','peakLagSec','peakCorr','chunkHalf', ...
                'doBaselineNorm','jobInd','permInd','shiftInd','jobType','baseDir', ...
                'runtimeSec','beh','nTrialsThisBehavior','nTrialsUsed','labelsToUse', ...
                'pyrAvgTrace','intAvgTrace','channelsToUseThisSess');
            continue;
        end

        %% ====================== PERM ======================
        shiftCorrUpper95 = getShiftCorrUpper95Behavior(outDir,labelTag,beh);
        if ~isfinite(shiftCorrUpper95)
            warning('No finite shift threshold for %s beh %d. Run shift jobs first.',baseDir,beh);
            continue;
        end

        acceptedPerm = false;
        tryCount = 0;
        xc = nan(1,numel(lags));
        peakLagSec = NaN;
        peakCorr = NaN;
        pyrAvgTrace = [];
        intAvgTrace = [];
        nTrialsThisBehavior = 0;
        nTrialsUsed = 0;

        tic;
        while ~acceptedPerm && tryCount < maxPermTries
            tryCount = tryCount + 1;

            [pyr_byCh,int_byCh,nTrials_try] = buildRealOrPermBehaviorSample( ...
                S,labelsPerTransition,beh,channelsToUseThisSess,tAxis,doBaselineNorm,true);
            if nTrials_try < 1, continue; end

            [pyrAvg_try,intAvg_try,nTrialsUsed_try] = averageAcrossTrialsAndChannels(pyr_byCh,int_byCh);
            if nTrialsUsed_try < 1, continue; end

            [xc_try,lag_try,corr_try] = lagSweepTrialAverage( ...
                pyrAvg_try,intAvg_try,lags,pyrStart,pyrEnd,T,binSize);

            if isfinite(corr_try) && isfinite(lag_try) && corr_try >= shiftCorrUpper95
                acceptedPerm = true;
                xc = xc_try;
                peakLagSec = lag_try;
                peakCorr = corr_try;
                pyrAvgTrace = pyrAvg_try;
                intAvgTrace = intAvg_try;
                nTrialsThisBehavior = nTrials_try;
                nTrialsUsed = nTrialsUsed_try;
            end
        end
        runtimeSec = toc;

        outFile = fullfile(outDir,sprintf( ...
            'concatCrossCorrPerCanonicalBehavior_%s_perm_%03d_beh%02d.mat',labelTag,permInd,beh));

        save(outFile,'lags','binSize','xc','peakLagSec','peakCorr','shiftCorrUpper95', ...
            'acceptedPerm','tryCount','maxPermTries','chunkHalf','doBaselineNorm', ...
            'jobInd','permInd','shiftInd','jobType','baseDir','runtimeSec','beh', ...
            'nTrialsThisBehavior','nTrialsUsed','labelsToUse','pyrAvgTrace', ...
            'intAvgTrace','channelsToUseThisSess');

        fprintf('    perm %d beh %d thresh=%.6f accepted=%d tries=%d corr=%.6f lag=%.6f\n', ...
            permInd,beh,shiftCorrUpper95,acceptedPerm,tryCount,peakCorr,peakLagSec);
    end
end

%% ================= helpers =================

function [pyr_byCh,int_byCh,nTrialsThisBehavior] = buildRealOrPermBehaviorSample( ...
    S,labelsPerTransition,beh,channelsToUseThisSess,tAxis,doBaselineNorm,doPermutation)

pyr_byCh = cell(1,numel(channelsToUseThisSess));
int_byCh = cell(1,numel(channelsToUseThisSess));
nTrialsThisBehavior = 0;

for ci = 1:numel(channelsToUseThisSess)
    ch = channelsToUseThisSess(ci);
    if ch > numel(labelsPerTransition) || isempty(labelsPerTransition{ch}), continue; end

    pyrWin = S.pyrCxWinCell{ch};
    intWin = S.intCxWinCell{ch};
    if isempty(pyrWin) || isempty(intWin), continue; end

    labels = labelsPerTransition{ch}(:);
    nEvt = min([size(pyrWin,1),size(intWin,1),numel(labels)]);
    if nEvt < 1, continue; end

    labels = labels(1:nEvt);
    evIdx = find(labels == beh);
    if isempty(evIdx), continue; end

    pyrWin_sel = pyrWin(evIdx,:,:);
    intWin_sel = intWin(evIdx,:,:);

    if doPermutation
        [nEvtP2,nPyr,~] = size(pyrWin_sel);
        [nEvtI2,nInt,~] = size(intWin_sel);
        nEvtUse = min(nEvtP2,nEvtI2);
        if nPyr+nInt < 2 || nEvtUse < 1, continue; end

        pooledWin = cat(2,pyrWin_sel(1:nEvtUse,:,:),intWin_sel(1:nEvtUse,:,:));
        idx = randperm(nPyr+nInt);
        idxP = idx(1:nPyr);
        idxI = idx(nPyr+1:end);
        pyrEvt = squeeze(mean(pooledWin(:,idxP,:),2,'omitnan'));
        intEvt = squeeze(mean(pooledWin(:,idxI,:),2,'omitnan'));
        if isvector(pyrEvt), pyrEvt = pyrEvt(:)'; end
        if isvector(intEvt), intEvt = intEvt(:)'; end
    else
        pyrEvt = meanEvt_keepEvents(pyrWin_sel);
        intEvt = meanEvt_keepEvents(intWin_sel);
    end

    if doBaselineNorm
        pyrEvt = subtractTrialBaseline(pyrEvt,tAxis,-500,-450);
        intEvt = subtractTrialBaseline(intEvt,tAxis,-500,-450);
    end

    pyr_byCh{ci} = pyrEvt;
    int_byCh{ci} = intEvt;
    nTrialsThisBehavior = nTrialsThisBehavior + size(pyrEvt,1);
end
end

function [pyr_byCh,int_byCh,nTrialsThisBehavior] = buildShiftedBehaviorSample( ...
    S,labelsPerTransition,beh,shiftCol,cortexIntFR,tAxis,doBaselineNorm, ...
    useFallbackShift,randShift,frameEMGSamples,frameNeuropixelSamples,totalLength,baseDir,channelsToUseThisSess)

pyr_byCh = cell(1,numel(channelsToUseThisSess));
int_byCh = cell(1,numel(channelsToUseThisSess));
nTrialsThisBehavior = 0;

for ci = 1:numel(channelsToUseThisSess)
    ch = channelsToUseThisSess(ci);
    if ch > numel(labelsPerTransition) || isempty(labelsPerTransition{ch}), continue; end

    pyrWin = S.pyrCxWinCell{ch};
    intWin = S.intCxWinCell{ch};
    if isempty(pyrWin) || isempty(intWin), continue; end

    labels = labelsPerTransition{ch}(:);
    nEvt = min([size(pyrWin,1),size(intWin,1),numel(labels)]);
    if nEvt < 1, continue; end

    labels = labels(1:nEvt);
    evIdx = find(labels == beh);
    if isempty(evIdx), continue; end
    pyrWin_sel = pyrWin(evIdx,:,:);

    if ~useFallbackShift
        idx1k_all = S.validTransitionsNeurShiftedCell{ch,shiftCol};
        if isempty(idx1k_all), continue; end
        idx1k = idx1k_all(1:nEvt);
        idx1k = idx1k(evIdx);
    else
        validTransitions = S.validTransitionsCell{ch};
        if isempty(validTransitions), continue; end
        validTransitions = validTransitions(1:nEvt);
        shifted = validTransitions(evIdx) + randShift;
        shifted(shifted > totalLength) = shifted(shifted > totalLength) - (totalLength+1);
        shifted(shifted < 1) = shifted(shifted < 1) + totalLength;
        currentDir = pwd;
        cd(baseDir);
        neuralIndices = NeurEMGSync(shifted*20,frameEMGSamples,frameNeuropixelSamples,'EMG');
        cd(currentDir);
        idx1k = round(neuralIndices/30);
    end

    validEvtMask = ~isnan(idx1k);
    idx1k = idx1k(validEvtMask);
    pyrWin_sel = pyrWin_sel(validEvtMask,:,:);
    if isempty(idx1k) || isempty(pyrWin_sel), continue; end

    pyrEvt = meanEvt_keepEvents(pyrWin_sel);
    intEvt = meanEvt_keepEvents(extractWindowsFromFR(cortexIntFR,idx1k,tAxis));
    nEvtUse = min(size(pyrEvt,1),size(intEvt,1));
    if nEvtUse < 1, continue; end
    pyrEvt = pyrEvt(1:nEvtUse,:);
    intEvt = intEvt(1:nEvtUse,:);

    if doBaselineNorm
        pyrEvt = subtractTrialBaseline(pyrEvt,tAxis,-500,-450);
        intEvt = subtractTrialBaseline(intEvt,tAxis,-500,-450);
    end

    pyr_byCh{ci} = pyrEvt;
    int_byCh{ci} = intEvt;
    nTrialsThisBehavior = nTrialsThisBehavior + size(pyrEvt,1);
end
end

function Eout = subtractTrialBaseline(Ein,tAxisSec,tStartSec,tEndSec)
[~,iStart] = min(abs(tAxisSec-tStartSec));
[~,iEnd] = min(abs(tAxisSec-tEndSec));
if iEnd < iStart, [iStart,iEnd] = deal(iEnd,iStart); end
baselines = mean(Ein(:,iStart:iEnd),2,'omitnan');
Eout = Ein - baselines;
end

function E = meanEvt_keepEvents(M)
if isempty(M), E = []; return; end
E = mean(M,2,'omitnan');
E = permute(E,[1 3 2]);
E = reshape(E,size(E,1),size(E,2));
end

function [pyrAvgTrace,intAvgTrace,nTrialsUsed] = averageAcrossTrialsAndChannels(pyr_byCh,int_byCh)
pyrAll = [];
intAll = [];
for ci = 1:numel(pyr_byCh)
    pEvt = pyr_byCh{ci};
    iEvt = int_byCh{ci};
    if isempty(pEvt) || isempty(iEvt), continue; end
    nEvt = min(size(pEvt,1),size(iEvt,1));
    if nEvt < 1, continue; end
    pyrAll = cat(1,pyrAll,pEvt(1:nEvt,:));
    intAll = cat(1,intAll,iEvt(1:nEvt,:));
end
nTrialsUsed = size(pyrAll,1);
if nTrialsUsed < 1, pyrAvgTrace=[]; intAvgTrace=[]; return; end
pyrAvgTrace = mean(pyrAll,1,'omitnan');
intAvgTrace = mean(intAll,1,'omitnan');
end

function [xc,peakLagSec,peakCorr] = lagSweepTrialAverage(pyrAvgTrace,intAvgTrace,lags,pyrStart,pyrEnd,T,binSize)
xc = nan(1,numel(lags));
for iL = 1:numel(lags)
    L = lags(iL);
    intStart = pyrStart - L;
    intEnd = pyrEnd - L;
    if intStart < 1 || intEnd > T, continue; end
    pseg = pyrAvgTrace(pyrStart:pyrEnd);
    iseg = intAvgTrace(intStart:intEnd);
    valid = ~isnan(pseg) & ~isnan(iseg);
    if nnz(valid) > 10, xc(iL) = corr(pseg(valid)',iseg(valid)'); end
end
[peakCorr,peakIdx] = max(xc);
if isempty(peakIdx) || ~isfinite(peakCorr), peakLagSec=NaN; peakCorr=NaN; else, peakLagSec=lags(peakIdx)*binSize; end
end

function xcZeroLag = zeroLagTrialAverage(pyrAvgTrace,intAvgTrace,pyrStart,pyrEnd)
pseg = pyrAvgTrace(pyrStart:pyrEnd);
iseg = intAvgTrace(pyrStart:pyrEnd);
valid = ~isnan(pseg) & ~isnan(iseg);
xcZeroLag = NaN;
if nnz(valid) > 10, xcZeroLag = corr(pseg(valid)',iseg(valid)'); end
end

function shiftCorrUpper95 = getShiftCorrUpper95Behavior(outDir,labelTag,beh)
pattern = sprintf('concatCrossCorrPerCanonicalBehavior_%s_shift_*_zerolag_beh%02d.mat',labelTag,beh);
shiftFiles = dir(fullfile(outDir,pattern));
if isempty(shiftFiles), shiftCorrUpper95=NaN; return; end
shiftVals = nan(numel(shiftFiles),1);
for k = 1:numel(shiftFiles)
    D = load(fullfile(outDir,shiftFiles(k).name),'xcZeroLag');
    if isfield(D,'xcZeroLag') && isfinite(D.xcZeroLag), shiftVals(k)=D.xcZeroLag; end
end
shiftVals = shiftVals(isfinite(shiftVals));
if isempty(shiftVals), shiftCorrUpper95=NaN; return; end
shiftCorrUpper95 = prctile(shiftVals,95);
fprintf('    beh %d shiftCorrUpper95=%.6f from %d controls\n',beh,shiftCorrUpper95,numel(shiftVals));
end

function W = extractWindowsFromFR(frMat,neuralIdx1k,tAxis)
nEvt=numel(neuralIdx1k); nNeur=size(frMat,1); nTime=numel(tAxis); Ttot=size(frMat,2);
W=nan(nEvt,nNeur,nTime);
for e=1:nEvt
    t0=neuralIdx1k(e);
    if isnan(t0), continue; end
    rngIdx=t0+tAxis;
    if any(rngIdx<1)||any(rngIdx>Ttot), continue; end
    W(e,:,:)=frMat(:,rngIdx);
end
end

function [cortexIntFR,ok] = loadCortexIntFR(folderPath)
ok=false; cortexIntFR=[];
try
    frFile=fullfile(folderPath,'NeuralFiringRates1msBins10msGauss.mat');
    if ~isfile(frFile), return; end
    load(frFile,'cortexFRs','cortexInds');
    clsFile=fullfile('/home/asa7288/Transfer','AA_classifications.mat');
    if ~isfile(clsFile), return; end
    load(clsFile,'classifications');
    if contains(folderPath,'D026'), matchRow=1;
    elseif contains(folderPath,'D020'), matchRow=2;
    elseif contains(folderPath,'D024'), matchRow=3;
    elseif contains(folderPath,'D043'), matchRow=4;
    elseif contains(folderPath,'D050'), matchRow=5;
    elseif contains(folderPath,'D054'), matchRow=6;
    else, return; end
    clsRow=classifications(matchRow,:);
    cortexIntFR=cortexFRs(clsRow{1,1}(cortexInds)==1,:);
    ok=true;
catch ME
    warning('loadCortexIntFR failed for %s: %s',folderPath,ME.message);
    ok=false;
end
end
