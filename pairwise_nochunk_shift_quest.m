% pairwise_nochunk_shift_quest.m
% QUEST version: shifted no-chunk pairwise xcorr
% ONE JOB = ONE SHIFT (jobInd = 1..100)

clearvars; clc;

%% ---------- settings ----------
sessInd = 1;      % 1=D026, 2=D020, 3=D024, 4=D043
jobInd  = 1;      % shift index (1..100)

baseDirs = {
    '/home/asa7288/Transfer/D026', ...
    '/home/asa7288/Transfer/D020', ...
    '/home/asa7288/Transfer/D024', ...
    '/home/asa7288/Transfer/D043'
};

baseDir = baseDirs{sessInd};

%% ---------- parameters ----------
binSize = 0.001;

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

fprintf('sessInd=%d | shift jobInd=%d | nInt=%d | nPyr=%d\n', sessInd, jobInd, numInter, numPyr);

%% ---------- apply circular shift to ALL interneurons ----------
numBins = size(interFRs,2);
minShift = round(30/binSize);
maxShift = numBins - minShift;

rng(jobInd);  % reproducible per shift
shiftAmt = randi([minShift maxShift]);

interFRs = circshift(interFRs,[0 shiftAmt]);

%% ---------- compute 0-lag only ----------
nullCorrMat = nan(numInter, numPyr);

tic
for i = 1:numInter
    for j = 1:numPyr
        intTS = interFRs(i,:);
        pyrTS = pyrFRs(j,:);

        valid = ~isnan(intTS) & ~isnan(pyrTS);
        if nnz(valid) > 2
            nullCorrMat(i,j) = corr(intTS(valid)', pyrTS(valid)');
        end
    end
end
fprintf('shift %d done in %.1f s\n', jobInd, toc);

%% ---------- save ----------
outDir = fullfile(baseDir,'quest_runs');
if ~exist(outDir,'dir'), mkdir(outDir); end

outFile = fullfile(outDir, sprintf('pairwise_nochunk_sess%02d_shift_%03d.mat', sessInd, jobInd));

save(outFile, 'nullCorrMat', 'binSize', 'jobInd', 'sessInd', 'baseDir', 'shiftAmt', 'minShift', 'maxShift', 'numInter', 'numPyr', 'numBins', '-v7.3');

fprintf('saved %s\n', outFile);
