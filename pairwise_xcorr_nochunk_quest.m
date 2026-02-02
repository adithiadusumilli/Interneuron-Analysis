% pairwise_xcorr_nochunk_quest.m
% quest version of no-chunk pairwise cross-correlation
% one job = one session × one shift condition

% behavior:
%   - jobInd = 0  : full lag sweep for every (int,pyr) pair (xc curve + peak lag/corr)
%   - jobInd > 0  : shifted control, NO lag sweep (0-lag corr only per pair)
%                  (stores 0-lag in xcMat(:,:,1); other lag bins remain NaN)

%% ---------- job settings ----------
isShifted = jobInd > 0;

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

%% ---------- load classifications ----------
C = load('/home/asa7288/Transfer/AA_classifications.mat','classifications');
classifications = C.classifications;

if contains(baseDir,'D026'), matchRow = 1;
elseif contains(baseDir,'D020'), matchRow = 2;
elseif contains(baseDir,'D024'), matchRow = 3;
elseif contains(baseDir,'D043'), matchRow = 4;
else, error('unknown session folder'); end

%% ---------- load FR data ----------
F = load(fullfile(baseDir,'NeuralFiringRates1msBins10msGauss.mat'), 'cortexFRs','cortexInds');

frMatrix  = F.cortexFRs;
regionInds = F.cortexInds;

regionClass = classifications{matchRow,1}(regionInds);

interFRs = frMatrix(regionClass==1,:);
pyrFRs   = frMatrix(regionClass==0,:);

numInter = size(interFRs,1);
numPyr   = size(pyrFRs,1);

if numInter==0 || numPyr==0
    error('no neurons after classification');
end

%% ---------- shifted control ----------
if isShifted
    numBins = size(interFRs,2);
    minShift = round(30/binSize);
    maxShift = numBins - minShift;

    rng(jobInd); % reproducible per shift job
    shiftAmt = randi([minShift maxShift]);

    interFRs = circshift(interFRs,[0 shiftAmt]);
end

%% ---------- pairwise xcorr ----------
nL = numel(lags);

xcMat = nan(numInter,numPyr,nL);
peakLagMat = nan(numInter,numPyr);
peakCorrMat = nan(numInter,numPyr);

tic

if ~isShifted
    % ===========================
    % jobInd = 0 : full lag sweep
    % ===========================
    for i = 1:numInter
        intTS = interFRs(i,:);

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
                    xc(li) = corr(intSeg(valid)',pyrSeg(valid)');
                end
            end

            xcMat(i,j,:) = xc;
            [pk,idx] = max(xc);
            peakCorrMat(i,j) = pk;
            peakLagMat(i,j)  = lags(idx)*binSize;
        end
    end

else
    % ===========================
    % jobInd > 0 : 0-lag only
    % ===========================
    for i = 1:numInter
        intTS = interFRs(i,:);

        for j = 1:numPyr
            pyrTS = pyrFRs(j,:);

            valid = ~isnan(intTS) & ~isnan(pyrTS);
            if nnz(valid) > 2
                r = corr(intTS(valid)', pyrTS(valid)');
                xcMat(i,j,1)     = r;   % store as "0-lag" slice
                peakCorrMat(i,j) = r;
                peakLagMat(i,j)  = 0;   % by definition for 0-lag-only controls
            end
        end
    end
end

fprintf('pairwise done in %.1f s\n', toc);

%% ---------- save ----------
outDir = fullfile(baseDir,'quest_runs_nochunk');
if ~exist(outDir,'dir'), mkdir(outDir); end

if ~isShifted
    outFile = fullfile(outDir,'pairwise_nochunk_real.mat');
else
    outFile = fullfile(outDir, sprintf('pairwise_nochunk_shift_%03d.mat',jobInd));
end

save(outFile, 'xcMat','peakLagMat','peakCorrMat','lags','binSize', 'jobInd','baseDir','-v7.3');

fprintf('saved %s\n', outFile);
