% pairwise_nochunk_allPairs_real_or_shift_quest.m

% no-chunk ALL-PAIRS pairwise xcorr
% NO EMG windows, NO trial averaging

% jobInd = 0 -> real all-pairs full lag sweep
% jobInd = 1:100 -> shifted null matrix

% real:
%   computes xcorr for all i<j pairs

% shift:
%   for each i<j, corr(unshifted i vs shifted j), then mirror

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
    'define jobInd before running. e.g., matlab -batch "jobInd=0; pairwise_nochunk_allPairs_real_or_shift_quest"');

assert(jobInd >= 0 && jobInd <= 100, ...
    'jobInd must be 0 for real or 1:100 for shift null');

isShiftJob = jobInd > 0;

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
pyrFRs = frMatrix(regionClass == 0,:);

nInt = size(interFRs,1);
nPyr = size(pyrFRs,1);
nAll = nInt + nPyr;

typeVec = [ones(nInt,1); zeros(nPyr,1)]; % 1=int, 0=pyr
allFRs = [interFRs; pyrFRs]; % nAll x time
numBins = size(allFRs,2);

outDir = fullfile(baseDir,'quest_runs');
if ~exist(outDir,'dir')
    mkdir(outDir);
end

fprintf('\n=== NO-CHUNK ALL-PAIRS | sessInd=%d | jobInd=%d ===\n', sessInd, jobInd);
fprintf('baseDir=%s\n', baseDir);
fprintf('nInt=%d | nPyr=%d | nAll=%d | nLags=%d | numBins=%d\n', nInt, nPyr, nAll, nL, numBins);

%% ============================================================
% jobInd = 0: real all-pairs full xcorr
%% ============================================================
if ~isShiftJob

    xcMatAll = nan(nAll, nAll, nL);
    peakCorrMatAll = nan(nAll, nAll);
    peakLagMatAll = nan(nAll, nAll);

    tic;

    for i = 1:nAll
        ts_i = allFRs(i,:);

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

            xcMatAll(i,j,:) = xc;

            [pk, idx] = max(xc);
            peakCorrMatAll(i,j) = pk;
            peakLagMatAll(i,j) = lags(idx) * binSize;
        end
    end

    runtimeSec = toc;

    outFile = fullfile(outDir, ...
        sprintf('pairwise_nochunk_allPairs_sess%02d_real_fullxc.mat', sessInd));

    save(outFile, 'xcMatAll', 'peakCorrMatAll', 'peakLagMatAll', 'lags', 'binSize', ...
        'jobInd', 'sessInd', 'baseDir', 'matchRow', 'nInt', 'nPyr', 'nAll', 'typeVec', ...
        'numBins', 'runtimeSec', '-v7.3');

    fprintf('saved real all-pairs file:\n%s\n', outFile);
    fprintf('runtime %.1f s\n', runtimeSec);

%% ============================================================
% jobInd = 1:100: shifted null matrix
%% ============================================================
else

    rng(jobInd);

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
        sprintf('pairwise_nochunk_allPairs_sess%02d_shift_%03d.mat', sessInd, jobInd));

    save(outFile, 'nullCorrMatAll', 'binSize', 'jobInd', 'sessInd', 'baseDir', ...
        'matchRow', 'nInt', 'nPyr', 'nAll', 'typeVec', 'numBins', 'shiftAmtPerNeuron', 'minShiftBins', ...
        'runtimeSec', '-v7.3');

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
