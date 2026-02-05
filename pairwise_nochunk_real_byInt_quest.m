% pairwise_nochunk_real_byInt_quest.m
% QUEST version: unshifted no-chunk pairwise xcorr
% ONE JOB = ONE INTERNEURON -- parallel across interneurons yay
% --> jobInd IS the interneuron index

clearvars; clc;

%% ---------- settings ----------
sessInd = 1;      % 1=D026, 2=D020, 3=D024, 4=D043
jobInd  = 1;      % <-- THIS IS NOW INT INDEX !!!!!!!!!!!!

baseDirs = {
    '/home/asa7288/Transfer/D026', ...
    '/home/asa7288/Transfer/D020', ...
    '/home/asa7288/Transfer/D024', ...
    '/home/asa7288/Transfer/D043'
};

baseDir = baseDirs{sessInd};

%% ---------- parameters ----------
binSize = 0.001;
maxLagSecs = 0.5;
maxLagBins = round(maxLagSecs/binSize);
lags = -maxLagBins:maxLagBins;
nL = numel(lags);

%% ---------- load classifications ----------
C = load('/home/asa7288/Transfer/AA_classifications.mat','classifications');
classifications = C.classifications;

if contains(baseDir,'D026'), matchRow = 1;
elseif contains(baseDir,'D020'), matchRow = 2;
elseif contains(baseDir,'D024'), matchRow = 3;
elseif contains(baseDir,'D043'), matchRow = 4;
else, error('unknown session folder'); end

%% ---------- load FR data ----------
F = load(fullfile(baseDir,'NeuralFiringRates1msBins10msGauss.mat'), ...
    'cortexFRs','cortexInds');

frMatrix   = F.cortexFRs;
regionInds = F.cortexInds;

regionClass = classifications{matchRow,1}(regionInds);

interFRs = frMatrix(regionClass==1,:);
pyrFRs   = frMatrix(regionClass==0,:);

numInter = size(interFRs,1);
numPyr   = size(pyrFRs,1);

assert(jobInd >= 1 && jobInd <= numInter, ...
    'jobInd must be between 1 and numInter');

fprintf('sessInd=%d | baseDir=%s\n', sessInd, baseDir);
fprintf('jobInd (interneuron) = %d / %d | nPyr=%d | nLags=%d\n', ...
    jobInd, numInter, numPyr, nL);

%% ---------- compute for THIS interneuron only ----------
intTS = interFRs(jobInd,:);

xcRow      = nan(numPyr, nL);   % (pyr x lags)
peakLagRow = nan(numPyr, 1);
peakCorrRow = nan(numPyr, 1);

tic
for j = 1:numPyr
    pyrTS = pyrFRs(j,:);
    xc = nan(1,nL);

    for li = 1:nL
        L = lags(li);

        if L < 0
            intSeg = intTS(1:end+L);
            pyrSeg = pyrTS(1-L:end);
        elseif L > 0
            intSeg = intTS(1+L:end);
            pyrSeg = pyrTS(1:end-L);
        else
            intSeg = intTS;
            pyrSeg = pyrTS;
        end

        valid = ~isnan(intSeg) & ~isnan(pyrSeg);
        if nnz(valid) > 2
            xc(li) = corr(intSeg(valid)', pyrSeg(valid)');
        end
    end

    xcRow(j,:) = xc;
    [pk, idx] = max(xc);
    peakCorrRow(j) = pk;
    peakLagRow(j) = lags(idx)*binSize;
end
fprintf('interneuron %d done in %.1f s\n', jobInd, toc);

%% ---------- save ----------
outDir = fullfile(baseDir,'quest_runs_nochunk');
if ~exist(outDir,'dir'), mkdir(outDir); end

outFile = fullfile(outDir, sprintf('pairwise_nochunk_sess%02d_real_int%03d.mat', sessInd, jobInd));

runtimeSec = toc;  % moving toc into a var instead of just fprintf
timestamp = datestr(now,'yyyymmdd_HHMMSS');
nCortexCells = numel(regionInds);

save(outFile, 'xcRow','peakLagRow','peakCorrRow', 'lags','binSize', 'jobInd','sessInd','baseDir', 'numInter','numPyr', 'matchRow','nCortexCells', 'runtimeSec', 'timestamp', '-v7.3');

fprintf('saved %s\n', outFile);
