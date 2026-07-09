% pairwise_nochunk_allPairs_real_or_shift_quest.m
% no-chunk ALL-PAIRS pairwise xcorr
% NO EMG windows, NO trial averaging
%
% jobInd = 1:nAll  -> real row job; computes neuron i vs j>i
% jobInd = 101:200 -> shifted null matrix; shiftInd = jobInd - 100

clc;

%% ---------- settings ----------

sessInd = 4;  % 1=D026,2=D020,3=D024,4=D043,5=D050,6=D054

baseDirs = {
    '/home/asa7288/Transfer/D026', ...
    '/home/asa7288/Transfer/D020', ...
    '/home/asa7288/Transfer/D024', ...
    '/home/asa7288/Transfer/D043', ...
    '/home/asa7288/Transfer/D050', ...
    '/home/asa7288/Transfer/D054'
};

baseDir = baseDirs{sessInd};

%% ---------- quest: require jobInd ----------
assert(exist('jobInd','var')==1, ...
    'define jobInd before running. e.g., matlab -batch "jobInd=1; pairwise_nochunk_allPairs_real_or_shift_quest"');

%% ---------- parameters ----------
binSize = 0.001;
maxLagSecs = 0.5;
maxLagBins = round(maxLagSecs/binSize);
lags = -maxLagBins:maxLagBins;
nL = numel(lags);

minShiftSec = 0.03;  % >=30 ms
minShiftBins = round(minShiftSec/binSize);

%% ---------- load classifications ----------
C = load('/home/asa7288/Transfer/AA_classifications.mat','classifications');
classifications = C.classifications;

matchRow = getMatchRow(baseDir);

%% ---------- load FR data cortex only ----------
F = load(fullfile(baseDir,'NeuralFiringRates1msBins10msGauss.mat'), ...
    'cortexFRs','cortexInds');

frMatrix = F.cortexFRs;
regionInds = F.cortexInds;

regionClass = classifications{matchRow,1}(regionInds);

interFRs = frMatrix(regionClass == 1,:);
pyrFRs   = frMatrix(regionClass == 0,:);

nInt = size(interFRs,1);
nPyr = size(pyrFRs,1);
nAll = nInt + nPyr;

typeVec = [ones(nInt,1); zeros(nPyr,1)]; % 1=int, 0=pyr
allFRs = [interFRs; pyrFRs];             % nAll x time
numBins = size(allFRs,2);

outDir = fullfile(baseDir,'quest_runs');
if ~exist(outDir,'dir')
    mkdir(outDir);
end

%% ---------- decide job type after nAll is known ----------
if jobInd >= 1 && jobInd <= nAll
    jobType = "real";
    realRow = jobInd;
    shiftInd = NaN;
elseif jobInd >= 101 && jobInd <= 200
    jobType = "shift";
    realRow = NaN;
    shiftInd = jobInd - 100;
else
    error('jobInd must be 1:nAll for real rows or 101:200 for shift jobs. Here nAll=%d.', nAll);
end

fprintf('\n=== NO-CHUNK ALL-PAIRS | sessInd=%d | jobInd=%d | jobType=%s ===\n', ...
    sessInd, jobInd, jobType);
fprintf('baseDir=%s\n', baseDir);
fprintf('nInt=%d | nPyr=%d | nAll=%d | nLags=%d | numBins=%d\n', ...
    nInt, nPyr, nAll, nL, numBins);

%% ============================================================
% jobInd = 1:nAll: real row full xcorr
%% ============================================================
if jobType == "real"

    i = realRow;
    ts_i = allFRs(i,:);

    xcRowAll = nan(nAll, nL);        % only fills j>i
    peakCorrRowAll = nan(nAll, 1);
    peakLagRowAll = nan(nAll, 1);

    tic;

    for j = (i+1):nAll
        ts_j = allFRs(j,:);

        xc = nan(1,nL);

        for li = 1:nL
            L = lags(li);

            if L < 0
                a = ts_i(1:end+L);
                b = ts_j(1-L:end);
            elseif L > 0
                a = ts_i(1+L:end);
                b = ts_j(1:end-L);
            else
                a = ts_i;
                b = ts_j;
            end

            valid = ~isnan(a) & ~isnan(b);

            if nnz(valid) > 2
                xc(li) = corr(a(valid)', b(valid)');
            end
        end

        xcRowAll(j,:) = xc;

        [pk, idx] = max(xc);
        peakCorrRowAll(j) = pk;
        peakLagRowAll(j) = lags(idx) * binSize;
    end

    runtimeSec = toc;

    outFile = fullfile(outDir, ...
        sprintf('pairwise_nochunk_allPairs_sess%02d_real_row%03d.mat', sessInd, i));

    save(outFile, ...
        'xcRowAll', 'peakCorrRowAll', 'peakLagRowAll', ...
        'lags', 'binSize', ...
        'jobInd', 'realRow', 'shiftInd', 'jobType', ...
        'sessInd', 'baseDir', 'matchRow', ...
        'nInt', 'nPyr', 'nAll', 'typeVec', ...
        'numBins', 'runtimeSec', ...
        '-v7.3');

    fprintf('saved real row file:\n%s\n', outFile);
    fprintf('runtime %.1f s\n', runtimeSec);

%% ============================================================
% jobInd = 101:200: shifted null matrix
%% ============================================================
else

    rng(shiftInd);

    maxShiftBins = numBins - minShiftBins;
    assert(maxShiftBins > minShiftBins, 'not enough bins to shift safely');

    % random shift per neuron + random direction
    shiftAmtPerNeuron = randi([minShiftBins maxShiftBins], [nAll 1]);

    shiftSign = ones(nAll,1);
    shiftSign(rand(nAll,1) > 0.5) = -1;

    shiftAmtPerNeuron = shiftAmtPerNeuron .* shiftSign;

    allFRs_shift = allFRs;

    for k = 1:nAll
        allFRs_shift(k,:) = circshift(allFRs(k,:), [0 shiftAmtPerNeuron(k)]);
    end

    % null:
    % for each i<j, correlate unshifted i vs shifted j, then mirror
    nullCorrMatAll = nan(nAll, nAll);

    tic;

    for i = 1:nAll
        a = allFRs(i,:);  % unshifted i

        for j = (i+1):nAll
            b = allFRs_shift(j,:);  % shifted j

            valid = ~isnan(a) & ~isnan(b);

            if nnz(valid) > 2
                r = corr(a(valid)', b(valid)');
                nullCorrMatAll(i,j) = r;
                nullCorrMatAll(j,i) = r;
            end
        end
    end

    runtimeSec = toc;

    outFile = fullfile(outDir, ...
        sprintf('pairwise_nochunk_allPairs_sess%02d_shift_%03d.mat', sessInd, shiftInd));

    save(outFile, ...
        'nullCorrMatAll', ...
        'binSize', ...
        'jobInd', 'realRow', 'shiftInd', 'jobType', ...
        'sessInd', 'baseDir', 'matchRow', ...
        'nInt', 'nPyr', 'nAll', 'typeVec', ...
        'numBins', 'shiftAmtPerNeuron', 'minShiftBins', ...
        'runtimeSec', ...
        '-v7.3');

    fprintf('saved shifted null file:\n%s\n', outFile);
    fprintf('runtime %.1f s\n', runtimeSec);
end

%% ================= helper =================
function matchRow = getMatchRow(baseDir)

    if contains(baseDir,'D026')
        matchRow = 1;
    elseif contains(baseDir,'D020')
        matchRow = 2;
    elseif contains(baseDir,'D024')
        matchRow = 3;
    elseif contains(baseDir,'D043')
        matchRow = 4;
    elseif contains(baseDir,'D050')
        matchRow = 5;
    elseif contains(baseDir,'D054')
        matchRow = 6;
    else
        error('unknown session folder: %s', baseDir);
    end
end
