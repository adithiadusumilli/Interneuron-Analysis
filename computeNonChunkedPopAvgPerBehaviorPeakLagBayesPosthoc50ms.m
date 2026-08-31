function computeNonChunkedPopAvgPerBehaviorPeakLagBayesPosthoc50ms(savedCrossCorrFile)
% Post hoc PEAK-LAG model comparison for NON-CHUNKED population-average
% cross-correlation computed separately for each classifier behavior.

% Behavior-specific analogue of: computePopAvgPeakLagBayesPosthoc50ms

% Expected input: output from runNonChunkedXcorrPerBehavior_saveOnly_getMouseDataNames
%
% Expected structure:
%   results.sessions(i).mouseID
%   results.sessions(i).beh(b).beh
%   results.sessions(i).beh(b).peakLagSec
%   results.sessions(i).beh(b).permPeakLags

% If a behavior struct contains permAccepted, it is honored. Otherwise, finite permutation peak lags are treated as valid.

% Permutation blocks are preserved EXACTLY by original index:
%   1:100   = H0
%   101:200 = H+50
%   201:300 = H-50

% Post hoc transformations:
%   H0 = raw permutations 1:100
%   H+50 = raw permutations 101:200 + 0.050 s
%   H-50 = raw permutations 201:300 - 0.050 s

% H0 versus H+50 uses LEFT-tail empirical probabilities.
% H0 versus H-50 uses RIGHT-tail empirical probabilities.
% A +1 correction is applied to numerator and denominator.

% Evidence ratios:
%   H0/H+50 = pH0_left / pH50
%   H0/H-50 = pH0_right / pHneg50

% No sign flipping is applied.

% Default input: X:\David\AnalysesData\nonchunked_xcorr_by_classifier_cortex_all6Sessions_saved.mat

% Example: computeNonChunkedPopAvgPerBehaviorPeakLagBayesPosthoc50ms

arguments
    savedCrossCorrFile (1,1) string = ...
        "X:\David\AnalysesData\nonchunked_xcorr_by_classifier_cortex_all6Sessions_saved.mat"
end

%% ---------------- settings ----------------

lagShiftSec = 0.050;
makePlots = true;
nHistogramBins = 20;

permIndsH0 = 1:100;
permIndsH50 = 101:200;
permIndsHneg50 = 201:300;

%% ---------------- validate and load ----------------

if ~isfile(savedCrossCorrFile)
    error('Input file not found:\n%s',savedCrossCorrFile);
end

S = load(savedCrossCorrFile,'results');

if ~isfield(S,'results')
    error('Input file does not contain variable "results".');
end

results = S.results;

if ~isfield(results,'sessions') || isempty(results.sessions)
    error('results.sessions is missing or empty.');
end

sessions = results.sessions;
nSess = numel(sessions);

%% ---------------- determine behaviors ----------------

if isfield(results,'behaviors') && ~isempty(results.behaviors)
    behaviors = results.behaviors(:)';
else
    behaviors = [];

    for s = 1:nSess
        if isfield(sessions(s),'beh') && ~isempty(sessions(s).beh)
            behaviors = [sessions(s).beh.beh];
            break;
        end
    end
end

if isempty(behaviors)
    error('Could not determine behavior values from saved results.');
end

nBeh = numel(behaviors);

canonicalNames = { ...
    'unlabeled', ...
    'climbdown', ...
    'climbup', ...
    'eating', ...
    'grooming', ...
    'jumpdown', ...
    'jumping', ...
    'rearing', ...
    'still', ...
    'walkflat', ...
    'walkgrid'};

%% ---------------- initialize output ----------------

posthocResults = struct();
posthocResults.sourceFile = savedCrossCorrFile;
posthocResults.lagShiftSec = lagShiftSec;
posthocResults.metric = 'non-chunked population peak lag per behavior';
posthocResults.noSignFlipApplied = true;
posthocResults.behaviors = behaviors;

posthocResults.permIndsH0 = permIndsH0;
posthocResults.permIndsH50 = permIndsH50;
posthocResults.permIndsHneg50 = permIndsHneg50;

posthocResults.analysisDescription = [ ...
    'For each animal x behavior, the real non-chunked population-average ' ...
    'peak lag is compared with three independent permutation peak-lag ' ...
    'blocks. Permutations 1:100 are H0, permutations 101:200 are shifted ' ...
    'by +50 ms for H+50, and permutations 201:300 are shifted by -50 ms ' ...
    'for H-50. Original permutation indices are preserved and no sign ' ...
    'flipping is applied.'];

posthocResults.probabilityDefinition = [ ...
    'H0 versus H+50 uses left-tail empirical probabilities. H0 versus ' ...
    'H-50 uses right-tail empirical probabilities. A +1 correction is ' ...
    'applied to numerator and denominator.'];

emptyResult = struct( ...
    'mouseID','', ...
    'beh',NaN, ...
    'behaviorName','', ...
    'realPeakLagSec',NaN, ...
    'rawPeakLagsH0',[], ...
    'rawPeakLagsH50',[], ...
    'rawPeakLagsHneg50',[], ...
    'peakLagsH0',[], ...
    'peakLagsH50',[], ...
    'peakLagsHneg50',[], ...
    'nValidH0',0, ...
    'nValidH50',0, ...
    'nValidHneg50',0, ...
    'meanH0',NaN, ...
    'meanH50',NaN, ...
    'meanHneg50',NaN, ...
    'medianH0',NaN, ...
    'medianH50',NaN, ...
    'medianHneg50',NaN, ...
    'ci95H0',[NaN NaN], ...
    'ci95H50',[NaN NaN], ...
    'ci95Hneg50',[NaN NaN], ...
    'pValH0_left',NaN, ...
    'pValH0_right',NaN, ...
    'pValH50',NaN, ...
    'pValHneg50',NaN, ...
    'evidenceRatio_H0_over_H50',NaN, ...
    'evidenceRatio_H0_over_Hneg50',NaN);

posthocResults.sessions = repmat( ...
    struct('mouseID','','behaviors',[]),nSess,1);

%% ---------------- summary arrays ----------------

nRows = nSess*nBeh;

Animal = strings(nRows,1);
Behavior = nan(nRows,1);
BehaviorName = strings(nRows,1);
RealPeakLag_s = nan(nRows,1);

Nvalid_H0 = zeros(nRows,1);
Nvalid_H50 = zeros(nRows,1);
Nvalid_Hneg50 = zeros(nRows,1);

pValH0_left = nan(nRows,1);
pValH0_right = nan(nRows,1);
pValH50_left = nan(nRows,1);
pValHneg50_right = nan(nRows,1);

Evidence_H0_over_H50 = nan(nRows,1);
Evidence_H0_over_Hneg50 = nan(nRows,1);

row = 0;

%% ========================================================================
%% loop through animals
%% ========================================================================

for iSess = 1:nSess

    sess = sessions(iSess);

    if isfield(sess,'mouseID') && ~isempty(sess.mouseID)
        mouseID = string(sess.mouseID);
    else
        mouseID = "Animal_" + iSess;
    end

    posthocResults.sessions(iSess).mouseID = char(mouseID);
    behOut = repmat(emptyResult,nBeh,1);

    fprintf('\n============================================================\n');
    fprintf('ANIMAL: %s\n',mouseID);
    fprintf('============================================================\n');

    if ~isfield(sess,'beh') || isempty(sess.beh)
        warning('%s has no behavior results.',mouseID);
        posthocResults.sessions(iSess).behaviors = behOut;
        row = row + nBeh;
        continue;
    end

    sessBeh = sess.beh;

    %% ====================================================================
    %% loop through behaviors
    %% ====================================================================

    for bIdx = 1:nBeh

        row = row + 1;
        beh = behaviors(bIdx);

        Animal(row) = mouseID;
        Behavior(row) = beh;

        if beh >= 0 && beh <= 10
            behaviorName = string(canonicalNames{beh+1});
        else
            behaviorName = "behavior_" + beh;
        end

        BehaviorName(row) = behaviorName;

        R = emptyResult;
        R.mouseID = char(mouseID);
        R.beh = beh;
        R.behaviorName = char(behaviorName);

        fprintf('\n----------------------------------------\n');
        fprintf('%s | behavior %d (%s)\n',mouseID,beh,behaviorName);
        fprintf('----------------------------------------\n');

        thisIdx = find([sessBeh.beh] == beh,1,'first');

        if isempty(thisIdx)
            warning('%s has no result for behavior %d.',mouseID,beh);
            behOut(bIdx) = R;
            continue;
        end

        B = sessBeh(thisIdx);

        %% ---------- real peak lag ----------

        if ~isfield(B,'peakLagSec') || isempty(B.peakLagSec) || ...
                ~isscalar(B.peakLagSec) || ~isfinite(B.peakLagSec)

            warning('%s behavior %d has no valid peakLagSec.',mouseID,beh);
            behOut(bIdx) = R;
            continue;
        end

        realPeakLag = double(B.peakLagSec);
        R.realPeakLagSec = realPeakLag;
        RealPeakLag_s(row) = realPeakLag;

        %% ---------- 300 permutation slots ----------

        if ~isfield(B,'permPeakLags') || isempty(B.permPeakLags)

            warning('%s behavior %d is missing permPeakLags.',mouseID,beh);
            behOut(bIdx) = R;
            continue;
        end

        allRawPeakLags = double(B.permPeakLags(:)');

        if numel(allRawPeakLags) < 300
            warning(['%s behavior %d has only %d permutation peak-lag entries. ' ...
                     '300 are required for H0/H+50/H-50.'], ...
                     mouseID,beh,numel(allRawPeakLags));
            behOut(bIdx) = R;
            continue;
        end

        % Keep exactly the first 300 slots because their ORIGINAL index
        % defines model membership.
        allRawPeakLags = allRawPeakLags(1:300);

        %% ---------- accepted mask ----------

        if isfield(B,'permAccepted') && ~isempty(B.permAccepted) && ...
                numel(B.permAccepted) >= 300

            permAccepted = logical(B.permAccepted(:)');
            permAccepted = permAccepted(1:300);
        else
            % The saved non-chunked per-behavior structure may not contain
            % a separate accepted flag. In that case, finite lag entries
            % are the valid permutations.
            permAccepted = isfinite(allRawPeakLags);
        end

        %% ---------- extract independent blocks FIRST ----------

        rawH0 = allRawPeakLags(permIndsH0);
        rawH50 = allRawPeakLags(permIndsH50);
        rawHneg50 = allRawPeakLags(permIndsHneg50);

        acceptedH0 = permAccepted(permIndsH0);
        acceptedH50 = permAccepted(permIndsH50);
        acceptedHneg50 = permAccepted(permIndsHneg50);

        %% ---------- remove rejected / invalid within each block ----------

        rawPeakLagsH0 = rawH0(acceptedH0 & isfinite(rawH0));
        rawPeakLagsH50 = rawH50(acceptedH50 & isfinite(rawH50));
        rawPeakLagsHneg50 = rawHneg50(acceptedHneg50 & isfinite(rawHneg50));

        R.rawPeakLagsH0 = rawPeakLagsH0;
        R.rawPeakLagsH50 = rawPeakLagsH50;
        R.rawPeakLagsHneg50 = rawPeakLagsHneg50;

        if isempty(rawPeakLagsH0) || isempty(rawPeakLagsH50) || ...
                isempty(rawPeakLagsHneg50)

            warning('%s behavior %d has an empty valid model block.',mouseID,beh);
            behOut(bIdx) = R;
            continue;
        end

        %% ---------- construct model distributions ----------

        peakLagsH0 = rawPeakLagsH0;
        peakLagsH50 = rawPeakLagsH50 + lagShiftSec;
        peakLagsHneg50 = rawPeakLagsHneg50 - lagShiftSec;

        R.peakLagsH0 = peakLagsH0;
        R.peakLagsH50 = peakLagsH50;
        R.peakLagsHneg50 = peakLagsHneg50;

        R.nValidH0 = numel(peakLagsH0);
        R.nValidH50 = numel(peakLagsH50);
        R.nValidHneg50 = numel(peakLagsHneg50);

        Nvalid_H0(row) = R.nValidH0;
        Nvalid_H50(row) = R.nValidH50;
        Nvalid_Hneg50(row) = R.nValidHneg50;

        %% ---------- distribution summaries ----------

        [R.meanH0,R.medianH0,R.ci95H0] = ...
            summarizeDistribution(peakLagsH0);

        [R.meanH50,R.medianH50,R.ci95H50] = ...
            summarizeDistribution(peakLagsH50);

        [R.meanHneg50,R.medianHneg50,R.ci95Hneg50] = ...
            summarizeDistribution(peakLagsHneg50);

        %% ---------- empirical probabilities ----------

        R.pValH0_left = empiricalOneTailedProbability( ...
            realPeakLag,peakLagsH0,"left");

        R.pValH50 = empiricalOneTailedProbability( ...
            realPeakLag,peakLagsH50,"left");

        R.pValH0_right = empiricalOneTailedProbability( ...
            realPeakLag,peakLagsH0,"right");

        R.pValHneg50 = empiricalOneTailedProbability( ...
            realPeakLag,peakLagsHneg50,"right");

        %% ---------- evidence ratios ----------

        R.evidenceRatio_H0_over_H50 = ...
            safeRatio(R.pValH0_left,R.pValH50);

        R.evidenceRatio_H0_over_Hneg50 = ...
            safeRatio(R.pValH0_right,R.pValHneg50);

        pValH0_left(row) = R.pValH0_left;
        pValH0_right(row) = R.pValH0_right;
        pValH50_left(row) = R.pValH50;
        pValHneg50_right(row) = R.pValHneg50;

        Evidence_H0_over_H50(row) = R.evidenceRatio_H0_over_H50;
        Evidence_H0_over_Hneg50(row) = R.evidenceRatio_H0_over_Hneg50;

        %% ---------- print ----------

        fprintf('real peak lag: %+.3f s\n',realPeakLag);

        fprintf('valid model lags: H0=%d | H+50=%d | H-50=%d\n', ...
            R.nValidH0,R.nValidH50,R.nValidHneg50);

        fprintf('model medians: H0=%+.3f | H+50=%+.3f | H-50=%+.3f s\n', ...
            R.medianH0,R.medianH50,R.medianHneg50);

        fprintf('left-tail probabilities:  pH0=%.6f | pH+50=%.6f\n', ...
            R.pValH0_left,R.pValH50);

        fprintf('right-tail probabilities: pH0=%.6f | pH-50=%.6f\n', ...
            R.pValH0_right,R.pValHneg50);

        fprintf('evidence ratios: H0/H+50=%.6f | H0/H-50=%.6f\n', ...
            R.evidenceRatio_H0_over_H50, ...
            R.evidenceRatio_H0_over_Hneg50);

        %% ---------- plot ----------

        if makePlots
            plotAnimalBehaviorPeakLagDistributions( ...
                mouseID,beh,behaviorName,realPeakLag, ...
                peakLagsH0,peakLagsH50,peakLagsHneg50, ...
                nHistogramBins,R.pValH0_left,R.pValH0_right, ...
                R.pValH50,R.pValHneg50, ...
                R.evidenceRatio_H0_over_H50, ...
                R.evidenceRatio_H0_over_Hneg50);
        end

        behOut(bIdx) = R;
    end

    posthocResults.sessions(iSess).behaviors = behOut;
end

%% ---------------- summary table ----------------

summaryTable = table( ...
    Animal,Behavior,BehaviorName,RealPeakLag_s, ...
    Nvalid_H0,Nvalid_H50,Nvalid_Hneg50, ...
    pValH0_left,pValH0_right,pValH50_left,pValHneg50_right, ...
    Evidence_H0_over_H50,Evidence_H0_over_Hneg50, ...
    'VariableNames',{ ...
        'Animal', ...
        'Behavior', ...
        'BehaviorName', ...
        'RealPeakLag_s', ...
        'Nvalid_H0', ...
        'Nvalid_H50', ...
        'Nvalid_Hneg50', ...
        'pValH0_left', ...
        'pValH0_right', ...
        'pValH50_left', ...
        'pValHneg50_right', ...
        'Evidence_H0_over_H50', ...
        'Evidence_H0_over_Hneg50'});

posthocResults.summaryTable = summaryTable;

fprintf('\n============================================================\n');
fprintf('NON-CHUNKED POPAVG PER-BEHAVIOR BAYES SUMMARY\n');
fprintf('============================================================\n\n');

disp(summaryTable);

%% ---------------- save ----------------

[inputFolder,inputName,~] = fileparts(savedCrossCorrFile);

outFile = fullfile( ...
    inputFolder, ...
    inputName + "_PeakLagBayesPosthoc50ms.mat");

save(outFile,'posthocResults','summaryTable','-v7.3');

fprintf('\nsaved non-chunked per-behavior post hoc Bayes analysis to:\n%s\n', ...
    outFile);

end

%% ========================================================================
%% helper: empirical one-tailed probability
%% ========================================================================

function pVal = empiricalOneTailedProbability(realValue,modelDistribution,tail)

modelDistribution = modelDistribution(isfinite(modelDistribution));

if ~isfinite(realValue) || isempty(modelDistribution)
    pVal = NaN;
    return;
end

switch lower(string(tail))
    case "left"
        nExtreme = sum(modelDistribution <= realValue);
    case "right"
        nExtreme = sum(modelDistribution >= realValue);
    otherwise
        error('tail must be "left" or "right".');
end

pVal = (nExtreme + 1)/(numel(modelDistribution) + 1);

end

%% ========================================================================
%% helper: summarize distribution
%% ========================================================================

function [meanValue,medianValue,ci95] = summarizeDistribution(modelDistribution)

modelDistribution = modelDistribution(isfinite(modelDistribution));

if isempty(modelDistribution)
    meanValue = NaN;
    medianValue = NaN;
    ci95 = [NaN NaN];
    return;
end

meanValue = mean(modelDistribution,'omitnan');
medianValue = median(modelDistribution,'omitnan');

if numel(modelDistribution) >= 2
    ci95 = prctile(modelDistribution,[2.5 97.5]);
else
    ci95 = [NaN NaN];
end

end

%% ========================================================================
%% helper: safe ratio
%% ========================================================================

function ratio = safeRatio(numerator,denominator)

if ~isfinite(numerator) || ~isfinite(denominator) || denominator == 0
    ratio = NaN;
else
    ratio = numerator/denominator;
end

end

%% ========================================================================
%% helper: plot animal x behavior distributions
%% ========================================================================

function plotAnimalBehaviorPeakLagDistributions( ...
    mouseID,beh,behaviorName,realPeakLag, ...
    peakLagsH0,peakLagsH50,peakLagsHneg50, ...
    nHistogramBins,pValH0Left,pValH0Right,pValH50,pValHneg50, ...
    evidenceRatioH0H50,evidenceRatioH0Hneg50)

allValues = [ ...
    peakLagsH0(:); ...
    peakLagsH50(:); ...
    peakLagsHneg50(:); ...
    realPeakLag];

allValues = allValues(isfinite(allValues));

if isempty(allValues)
    return;
end

xMin = min(allValues);
xMax = max(allValues);

if xMin == xMax
    xPadding = 0.010;
else
    xPadding = 0.05*(xMax-xMin);
end

commonXLimits = [xMin-xPadding,xMax+xPadding];

commonEdges = linspace( ...
    commonXLimits(1),commonXLimits(2),nHistogramBins+1);

figure( ...
    'Name',sprintf('%s Beh %d Non-chunked Peak-Lag Models',mouseID,beh), ...
    'Color','w');

tl = tiledlayout(1,3, ...
    'TileSpacing','compact', ...
    'Padding','compact');

title(tl,sprintf( ...
    '%s | Behavior %d: %s | Non-chunked Population Peak-Lag Model Comparison', ...
    mouseID,beh,behaviorName), ...
    'FontSize',16);

%% H0
nexttile;
hold on;

histogram(peakLagsH0, ...
    'BinEdges',commonEdges, ...
    'FaceColor',[0.3 0.6 0.8], ...
    'EdgeColor','none');

xline(realPeakLag,'r-','LineWidth',2);
xlim(commonXLimits);

xlabel('Peak Lag (s)');
ylabel('Count');

title(sprintf('H0 | p_{left}=%.3f | p_{right}=%.3f', ...
    pValH0Left,pValH0Right));

box off;
set(gca,'FontSize',12,'TickDir','out');

%% H+50
nexttile;
hold on;

histogram(peakLagsH50, ...
    'BinEdges',commonEdges, ...
    'FaceColor',[0.3 0.6 0.8], ...
    'EdgeColor','none');

xline(realPeakLag,'r-','LineWidth',2);
xlim(commonXLimits);

xlabel('Peak Lag (s)');
ylabel('Count');

title(sprintf('H+50 | p_{left}=%.3f | H0/H+50=%.3f', ...
    pValH50,evidenceRatioH0H50));

box off;
set(gca,'FontSize',12,'TickDir','out');

%% H-50
nexttile;
hold on;

histogram(peakLagsHneg50, ...
    'BinEdges',commonEdges, ...
    'FaceColor',[0.3 0.6 0.8], ...
    'EdgeColor','none');

xline(realPeakLag,'r-','LineWidth',2);
xlim(commonXLimits);

xlabel('Peak Lag (s)');
ylabel('Count');

title(sprintf('H-50 | p_{right}=%.3f | H0/H-50=%.3f', ...
    pValHneg50,evidenceRatioH0Hneg50));

box off;
set(gca,'FontSize',12,'TickDir','out');

end
