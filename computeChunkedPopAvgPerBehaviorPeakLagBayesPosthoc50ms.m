function computeChunkedPopAvgPerBehaviorPeakLagBayesPosthoc50ms(savedCombinedFile)
% Post hoc PEAK-LAG model comparison for CHUNKED, TRIAL-AVERAGED, population-average cross-correlation computed separately by behavior.

% Behavior-specific analogue of computeChunkedPopAvgPeakLagBayesPosthoc50ms.
% Same logic:
%   perms 1:100 = H0
%   perms 101:200 = H+50, shifted +0.050 s post hoc
%   perms 201:300 = H-50, shifted -0.050 s post hoc
%   H0 vs H+50 uses LEFT tails
%   H0 vs H-50 uses RIGHT tails
%   empirical probabilities use +1 correction
%   evidence ratios are H0/H+50 and H0/H-50
%   no sign flipping

% Input must be the v2 combined behavior file with: allSessions.behaviorResults{bIdx} where each behavior result contains all animals.

if nargin < 1 || isempty(savedCombinedFile)
    savedCombinedFile = "C:\Users\mirilab\Documents\GlobusTransfer\concatCrossCorrPerCanonicalBehavior_classifier_chunked_trialavg_popavg_ALL_SESSIONS.mat";
end

lagShiftSec = 0.050;
makePlots = true;
nHistogramBins = 20;
permIndsH0 = 1:100;
permIndsH50 = 101:200;
permIndsHneg50 = 201:300;

if ~isfile(savedCombinedFile)
    error('Input file not found:\n%s', savedCombinedFile);
end

S = load(savedCombinedFile,'allSessions');
if ~isfield(S,'allSessions')
    error('Input file does not contain allSessions.');
end
A = S.allSessions;

if ~isfield(A,'behaviorResults') || isempty(A.behaviorResults)
    error('allSessions.behaviorResults missing/empty. Expected v2 combined structure.');
end
if ~isfield(A,'mouseIDs') || isempty(A.mouseIDs)
    error('allSessions.mouseIDs missing/empty.');
end

mouseIDs = string(A.mouseIDs(:));
nAnimals = numel(mouseIDs);

if isfield(A,'behaviors') && ~isempty(A.behaviors)
    behaviors = A.behaviors(:)';
else
    behaviors = nan(1,numel(A.behaviorResults));
    for bIdx = 1:numel(A.behaviorResults)
        behaviors(bIdx) = A.behaviorResults{bIdx}.beh;
    end
end
nBeh = numel(behaviors);

canonicalNames = {'climbdown','climbup','eating','grooming','jumpdown', ...
    'jumping','rearing','still','walkflat','walkgrid'};

posthocResults = struct();
posthocResults.sourceFile = savedCombinedFile;
posthocResults.metric = 'chunked trial-averaged population peak lag per behavior';
posthocResults.lagShiftSec = lagShiftSec;
posthocResults.noSignFlipApplied = true;
posthocResults.permIndsH0 = permIndsH0;
posthocResults.permIndsH50 = permIndsH50;
posthocResults.permIndsHneg50 = permIndsHneg50;
posthocResults.mouseIDs = mouseIDs;
posthocResults.behaviors = behaviors;
posthocResults.analysisDescription = ['For each animal x behavior, real peak lag is compared with three independent ' ...
    'random-label permutation peak-lag blocks. Permutations 1:100 are H0, 101:200 are shifted +50 ms, ' ...
    'and 201:300 are shifted -50 ms. Only accepted finite permutations are used, preserving original indices.'];

emptyResult = struct('mouseID','', 'beh',NaN, 'behaviorName','', ...
    'realPeakLagSec',NaN, 'rawPeakLagsH0',[], 'rawPeakLagsH50',[], 'rawPeakLagsHneg50',[], ...
    'peakLagsH0',[], 'peakLagsH50',[], 'peakLagsHneg50',[], ...
    'nValidH0',0, 'nValidH50',0, 'nValidHneg50',0, ...
    'nMissingOrRejectedH0',0, 'nMissingOrRejectedH50',0, 'nMissingOrRejectedHneg50',0, ...
    'meanH0',NaN, 'meanH50',NaN, 'meanHneg50',NaN, ...
    'medianH0',NaN, 'medianH50',NaN, 'medianHneg50',NaN, ...
    'ci95H0',[NaN NaN], 'ci95H50',[NaN NaN], 'ci95Hneg50',[NaN NaN], ...
    'pValH0_left',NaN, 'pValH0_right',NaN, 'pValH50',NaN, 'pValHneg50',NaN, ...
    'evidenceRatio_H0_over_H50',NaN, 'evidenceRatio_H0_over_Hneg50',NaN);

posthocResults.behaviorResults = cell(1,nBeh);

nRows = nAnimals*nBeh;
Animal = strings(nRows,1);
Behavior = nan(nRows,1);
BehaviorName = strings(nRows,1);
RealPeakLag_s = nan(nRows,1);
Nvalid_H0 = zeros(nRows,1);
Nvalid_H50 = zeros(nRows,1);
Nvalid_Hneg50 = zeros(nRows,1);
NmissingRejected_H0 = zeros(nRows,1);
NmissingRejected_H50 = zeros(nRows,1);
NmissingRejected_Hneg50 = zeros(nRows,1);
pValH0_left = nan(nRows,1);
pValH0_right = nan(nRows,1);
pValH50_left = nan(nRows,1);
pValHneg50_right = nan(nRows,1);
Evidence_H0_over_H50 = nan(nRows,1);
Evidence_H0_over_Hneg50 = nan(nRows,1);
row = 0;

for bIdx = 1:nBeh
    beh = behaviors(bIdx);
    B = A.behaviorResults{bIdx};
    if isempty(B), continue; end

    if beh >= 1 && beh <= numel(canonicalNames)
        behaviorName = string(canonicalNames{beh});
    else
        behaviorName = "behavior_" + beh;
    end

    fprintf('\n============================================================\n');
    fprintf('BEHAVIOR %d: %s\n',beh,behaviorName);
    fprintf('============================================================\n');

    behaviorOut = repmat(emptyResult,nAnimals,1);

    for a = 1:nAnimals
        row = row + 1;
        mouseID = mouseIDs(a);
        Animal(row) = mouseID;
        Behavior(row) = beh;
        BehaviorName(row) = behaviorName;

        R = emptyResult;
        R.mouseID = char(mouseID);
        R.beh = beh;
        R.behaviorName = char(behaviorName);

        if ~isfield(B,'real_peakLagSec') || numel(B.real_peakLagSec) < a || ~isfinite(B.real_peakLagSec(a))
            warning('%s behavior %d has no valid real peak lag.',mouseID,beh);
            behaviorOut(a) = R;
            continue;
        end

        realPeakLag = double(B.real_peakLagSec(a));
        R.realPeakLagSec = realPeakLag;
        RealPeakLag_s(row) = realPeakLag;

        if ~isfield(B,'perm_peakLagSec') || size(B.perm_peakLagSec,1) < a || size(B.perm_peakLagSec,2) < 300
            warning('%s behavior %d does not contain 300 permutation slots.',mouseID,beh);
            behaviorOut(a) = R;
            continue;
        end

        allRawPeakLags = double(B.perm_peakLagSec(a,1:300));

        if isfield(B,'perm_accepted') && size(B.perm_accepted,1) >= a && size(B.perm_accepted,2) >= 300
            acceptedMask = logical(B.perm_accepted(a,1:300));
        else
            acceptedMask = isfinite(allRawPeakLags);
            warning('%s behavior %d missing complete perm_accepted mask; using finite lags.',mouseID,beh);
        end
        validMask = acceptedMask & isfinite(allRawPeakLags);

        rawH0All = allRawPeakLags(permIndsH0);
        rawH50All = allRawPeakLags(permIndsH50);
        rawHneg50All = allRawPeakLags(permIndsHneg50);

        rawPeakLagsH0 = rawH0All(validMask(permIndsH0));
        rawPeakLagsH50 = rawH50All(validMask(permIndsH50));
        rawPeakLagsHneg50 = rawHneg50All(validMask(permIndsHneg50));

        R.rawPeakLagsH0 = rawPeakLagsH0;
        R.rawPeakLagsH50 = rawPeakLagsH50;
        R.rawPeakLagsHneg50 = rawPeakLagsHneg50;
        R.nValidH0 = numel(rawPeakLagsH0);
        R.nValidH50 = numel(rawPeakLagsH50);
        R.nValidHneg50 = numel(rawPeakLagsHneg50);
        R.nMissingOrRejectedH0 = 100 - R.nValidH0;
        R.nMissingOrRejectedH50 = 100 - R.nValidH50;
        R.nMissingOrRejectedHneg50 = 100 - R.nValidHneg50;

        Nvalid_H0(row)=R.nValidH0; Nvalid_H50(row)=R.nValidH50; Nvalid_Hneg50(row)=R.nValidHneg50;
        NmissingRejected_H0(row)=R.nMissingOrRejectedH0;
        NmissingRejected_H50(row)=R.nMissingOrRejectedH50;
        NmissingRejected_Hneg50(row)=R.nMissingOrRejectedHneg50;

        if isempty(rawPeakLagsH0) || isempty(rawPeakLagsH50) || isempty(rawPeakLagsHneg50)
            warning('%s behavior %d has an empty model block after filtering.',mouseID,beh);
            behaviorOut(a) = R;
            continue;
        end

        peakLagsH0 = rawPeakLagsH0;
        peakLagsH50 = rawPeakLagsH50 + lagShiftSec;
        peakLagsHneg50 = rawPeakLagsHneg50 - lagShiftSec;

        R.peakLagsH0 = peakLagsH0;
        R.peakLagsH50 = peakLagsH50;
        R.peakLagsHneg50 = peakLagsHneg50;

        [R.meanH0,R.medianH0,R.ci95H0] = summarizeDistribution(peakLagsH0);
        [R.meanH50,R.medianH50,R.ci95H50] = summarizeDistribution(peakLagsH50);
        [R.meanHneg50,R.medianHneg50,R.ci95Hneg50] = summarizeDistribution(peakLagsHneg50);

        R.pValH0_left = empiricalOneTailedProbability(realPeakLag,peakLagsH0,"left");
        R.pValH50 = empiricalOneTailedProbability(realPeakLag,peakLagsH50,"left");
        R.pValH0_right = empiricalOneTailedProbability(realPeakLag,peakLagsH0,"right");
        R.pValHneg50 = empiricalOneTailedProbability(realPeakLag,peakLagsHneg50,"right");

        R.evidenceRatio_H0_over_H50 = safeRatio(R.pValH0_left,R.pValH50);
        R.evidenceRatio_H0_over_Hneg50 = safeRatio(R.pValH0_right,R.pValHneg50);

        pValH0_left(row)=R.pValH0_left; pValH0_right(row)=R.pValH0_right;
        pValH50_left(row)=R.pValH50; pValHneg50_right(row)=R.pValHneg50;
        Evidence_H0_over_H50(row)=R.evidenceRatio_H0_over_H50;
        Evidence_H0_over_Hneg50(row)=R.evidenceRatio_H0_over_Hneg50;

        fprintf('%s | real=%+.3f s | valid H0/H+50/H-50 = %d/%d/%d | H0/H+50=%.4f | H0/H-50=%.4f\n', ...
            mouseID,realPeakLag,R.nValidH0,R.nValidH50,R.nValidHneg50, ...
            R.evidenceRatio_H0_over_H50,R.evidenceRatio_H0_over_Hneg50);

        if makePlots
            plotAnimalBehaviorPeakLagDistributions(mouseID,beh,behaviorName,realPeakLag, ...
                peakLagsH0,peakLagsH50,peakLagsHneg50,nHistogramBins, ...
                R.pValH0_left,R.pValH0_right,R.pValH50,R.pValHneg50, ...
                R.evidenceRatio_H0_over_H50,R.evidenceRatio_H0_over_Hneg50);
        end

        behaviorOut(a) = R;
    end

    posthocResults.behaviorResults{bIdx} = behaviorOut;
end

keep = 1:row;
summaryTable = table(Animal(keep),Behavior(keep),BehaviorName(keep),RealPeakLag_s(keep), ...
    Nvalid_H0(keep),Nvalid_H50(keep),Nvalid_Hneg50(keep), ...
    NmissingRejected_H0(keep),NmissingRejected_H50(keep),NmissingRejected_Hneg50(keep), ...
    pValH0_left(keep),pValH0_right(keep),pValH50_left(keep),pValHneg50_right(keep), ...
    Evidence_H0_over_H50(keep),Evidence_H0_over_Hneg50(keep), ...
    'VariableNames',{'Animal','Behavior','BehaviorName','RealPeakLag_s', ...
    'Nvalid_H0','Nvalid_H50','Nvalid_Hneg50', ...
    'NmissingRejected_H0','NmissingRejected_H50','NmissingRejected_Hneg50', ...
    'pValH0_left','pValH0_right','pValH50_left','pValHneg50_right', ...
    'Evidence_H0_over_H50','Evidence_H0_over_Hneg50'});

posthocResults.summaryTable = summaryTable;
fprintf('\n================ SUMMARY ================\n');
disp(summaryTable);

[inputFolder,inputName,~] = fileparts(savedCombinedFile);
outFile = fullfile(inputFolder,inputName + "_PeakLagBayesPosthoc50ms.mat");
save(outFile,'posthocResults','summaryTable','-v7.3');
fprintf('\nsaved behavior-specific post hoc Bayes analysis to:\n%s\n',outFile);
end

function pVal = empiricalOneTailedProbability(realValue,modelDistribution,tail)
modelDistribution = modelDistribution(isfinite(modelDistribution));
if ~isfinite(realValue) || isempty(modelDistribution), pVal=NaN; return; end
switch lower(string(tail))
    case "left", nExtreme = sum(modelDistribution <= realValue);
    case "right", nExtreme = sum(modelDistribution >= realValue);
    otherwise, error('tail must be "left" or "right".');
end
pVal = (nExtreme + 1)/(numel(modelDistribution) + 1);
end

function [meanValue,medianValue,ci95] = summarizeDistribution(x)
x = x(isfinite(x));
if isempty(x), meanValue=NaN; medianValue=NaN; ci95=[NaN NaN]; return; end
meanValue = mean(x,'omitnan');
medianValue = median(x,'omitnan');
if numel(x)>=2, ci95=prctile(x,[2.5 97.5]); else, ci95=[NaN NaN]; end
end

function ratio = safeRatio(numerator,denominator)
if ~isfinite(numerator) || ~isfinite(denominator) || denominator==0
    ratio=NaN;
else
    ratio=numerator/denominator;
end
end

function plotAnimalBehaviorPeakLagDistributions(mouseID,beh,behaviorName,realPeakLag, ...
    peakLagsH0,peakLagsH50,peakLagsHneg50,nHistogramBins, ...
    pValH0Left,pValH0Right,pValH50,pValHneg50,evidenceRatioH0H50,evidenceRatioH0Hneg50)

allValues=[peakLagsH0(:);peakLagsH50(:);peakLagsHneg50(:);realPeakLag];
allValues=allValues(isfinite(allValues));
if isempty(allValues), return; end
xMin=min(allValues); xMax=max(allValues);
if xMin==xMax, pad=0.010; else, pad=0.05*(xMax-xMin); end
xLim=[xMin-pad xMax+pad];
edges=linspace(xLim(1),xLim(2),nHistogramBins+1);

figure('Name',sprintf('%s Beh %d Peak-Lag Models',mouseID,beh),'Color','w');
tl=tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
title(tl,sprintf('%s | Behavior %d: %s | Chunked Popavg Peak-Lag Model Comparison', ...
    mouseID,beh,behaviorName),'FontSize',16);

nexttile; hold on;
histogram(peakLagsH0,'BinEdges',edges,'FaceColor',[0.3 0.6 0.8],'EdgeColor','none');
xline(realPeakLag,'r-','LineWidth',2); xlim(xLim); xlabel('Peak Lag (s)'); ylabel('Count');
title(sprintf('H0 | p_{left}=%.3f | p_{right}=%.3f',pValH0Left,pValH0Right)); box off;

nexttile; hold on;
histogram(peakLagsH50,'BinEdges',edges,'FaceColor',[0.3 0.6 0.8],'EdgeColor','none');
xline(realPeakLag,'r-','LineWidth',2); xlim(xLim); xlabel('Peak Lag (s)'); ylabel('Count');
title(sprintf('H+50 | p_{left}=%.3f | H0/H+50=%.3f',pValH50,evidenceRatioH0H50)); box off;

nexttile; hold on;
histogram(peakLagsHneg50,'BinEdges',edges,'FaceColor',[0.3 0.6 0.8],'EdgeColor','none');
xline(realPeakLag,'r-','LineWidth',2); xlim(xLim); xlabel('Peak Lag (s)'); ylabel('Count');
title(sprintf('H-50 | p_{right}=%.3f | H0/H-50=%.3f',pValHneg50,evidenceRatioH0Hneg50)); box off;
end
