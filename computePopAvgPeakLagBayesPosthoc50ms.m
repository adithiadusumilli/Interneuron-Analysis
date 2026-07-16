function computePopAvgPeakLagBayesPosthoc50ms(savedCrossCorrFile)
% post hoc population-averaged peak-lag model comparison
% This analysis uses peak lag itself as the model-comparison metric.

% Expected input:
%   output from: runCrossCorrelation_saveOnly_getMouseDataNames

% The input file must contain:
%   xcorrResults.sessions(i).peakLag
%   xcorrResults.sessions(i).permPeakLagsAll300
%   xcorrResults.sessions(i).permAccepted

% Permutation blocks:
%     1:100 = H0 raw peak lags
%   101:200 = H+50 raw peak lags
%   201:300 = H-50 raw peak lags

% Post hoc transformations:
%   H0 peak lags    = raw permutations 1:100
%   H+50 peak lags  = raw permutations 101:200 + 0.050 s
%   H-50 peak lags  = raw permutations 201:300 - 0.050 s

% For each model, the real peak lag is compared with the model-specific peak-lag distribution using an empirical two-sided compatibility probability centered on that distribution's median:
%   pModel =
%       fraction of null deviations from the model median that are
%       at least as large as the real lag's deviation from that median

% Larger pModel: real peak lag is more typical under that model

% Evidence ratios:
%   H0/H+50 = pH0 / pH50
%   H0/H-50 = pH0 / pHneg50

% Interpretation:
%   ratio > 1  favors H0
%   ratio < 1  favors the shifted model
%   ratio ~ 1  little distinction

% The function also creates one figure per animal containing the H0,
% H+50, and H-50 peak-lag frequency distributions.

% RUN:
% computePopAvgPeakLagBayesPosthoc50ms("X:\David\AnalysesData\InterneuronAnalyses\Lab Meeting Pres\6 animals run cross correlation, no chunking, pop-wise\runCrossCorrelation_savedOutputs_all6Animals.mat")

arguments
    savedCrossCorrFile (1,1) string
end

%% ---------------- settings ----------------

lagShiftSec = 0.050;
makePlots = true;
nHistogramBins = 20;

%% ---------------- validate and load ----------------

if ~isfile(savedCrossCorrFile)
    error('Input file not found:\n%s', savedCrossCorrFile);
end

S = load(savedCrossCorrFile, 'xcorrResults');

if ~isfield(S, 'xcorrResults')
    error('Input file does not contain xcorrResults.');
end

xcorrResults = S.xcorrResults;

if ~isfield(xcorrResults, 'sessions') || isempty(xcorrResults.sessions)
    error('xcorrResults.sessions is missing or empty.');
end

sessions = xcorrResults.sessions;
nSess = numel(sessions);

%% ---------------- determine permutation blocks ----------------

if isfield(xcorrResults, 'permIndsH0') && ~isempty(xcorrResults.permIndsH0)
    permIndsH0 = xcorrResults.permIndsH0;
else
    permIndsH0 = 1:100;
end

if isfield(xcorrResults, 'permIndsH50') && ~isempty(xcorrResults.permIndsH50)

    permIndsH50 = xcorrResults.permIndsH50;
else
    permIndsH50 = 101:200;
end

if isfield(xcorrResults, 'permIndsHneg50') && ~isempty(xcorrResults.permIndsHneg50)

    permIndsHneg50 = xcorrResults.permIndsHneg50;
else
    permIndsHneg50 = 201:300;
end

%% ---------------- initialize output ----------------

posthocResults = struct();
posthocResults.sourceFile = savedCrossCorrFile;
posthocResults.lagShiftSec = lagShiftSec;
posthocResults.metric = 'population peak lag';
posthocResults.analysisDescription = [ ...
    'The real population peak lag was compared with three independent ' ...
    'label-permutation peak-lag distributions. Permutations 1:100 were ' ...
    'used unshifted for H0, permutations 101:200 were shifted by +50 ms ' ...
    'for H+50, and permutations 201:300 were shifted by -50 ms for H-50.'];

posthocResults.probabilityDefinition = [ ...
    'Empirical two-sided compatibility probability based on the absolute ' ...
    'distance of the real peak lag from each model distribution median.'];

posthocResults.permIndsH0 = permIndsH0;
posthocResults.permIndsH50 = permIndsH50;
posthocResults.permIndsHneg50 = permIndsHneg50;

posthocResults.sessions = repmat(struct('mouseID', '', 'baseSessionName', '', 'processedDataFolder', '', ...
    'realPeakLagSec', NaN, 'realPeakCorr', NaN, 'rawPeakLagsH0', [], 'rawPeakLagsH50', [], 'rawPeakLagsHneg50', [], ...
    'peakLagsH0', [], 'peakLagsH50', [], 'peakLagsHneg50', [], 'nValidH0', 0, 'nValidH50', 0, 'nValidHneg50', 0, ...
    'medianH0', NaN, 'medianH50', NaN, 'medianHneg50', NaN, 'meanH0', NaN, 'meanH50', NaN, 'meanHneg50', NaN, 'ci95H0', [NaN NaN], ...
    'ci95H50', [NaN NaN], 'ci95Hneg50', [NaN NaN], 'pValH0', NaN, 'pValH50', NaN, 'pValHneg50', NaN, ...
    'evidenceRatio_H0_over_H50', NaN, 'evidenceRatio_H0_over_Hneg50', NaN), nSess, 1);

%% ---------------- summary arrays ----------------

animalIDs = strings(nSess,1);
realPeakLags = nan(nSess,1);
pH0All = nan(nSess,1);
pH50All = nan(nSess,1);
pHneg50All = nan(nSess,1);
evidenceH0overH50 = nan(nSess,1);
evidenceH0overHneg50 = nan(nSess,1);

%% ========================================================================
%% loop through animals
%% ========================================================================

for iSess = 1:nSess

    sess = sessions(iSess);

    %% ---------- animal ID ----------

    if isfield(sess, 'mouseID') && ~isempty(sess.mouseID)
        mouseID = string(sess.mouseID);
    elseif isfield(sess, 'animalLabel') && ~isempty(sess.animalLabel)
        mouseID = string(sess.animalLabel);
    else
        mouseID = "Animal_" + iSess;
    end

    animalIDs(iSess) = mouseID;

    fprintf('\n========================================\n');
    fprintf('processing %s\n', mouseID);
    fprintf('========================================\n');

    %% ---------- validate required fields ----------

    if ~isfield(sess, 'peakLag') || isempty(sess.peakLag) || ~isfinite(sess.peakLag)
        warning('%s has no valid real peak lag. Skipping.', mouseID);
        continue;
    end

    if ~isfield(sess, 'permPeakLagsAll300') || isempty(sess.permPeakLagsAll300)
        warning('%s is missing permPeakLagsAll300. Skipping.', mouseID);
        continue;
    end

    allRawPeakLags = sess.permPeakLagsAll300(:)';

    maxRequiredIndex = max([permIndsH0(:); permIndsH50(:); permIndsHneg50(:)]);

    if numel(allRawPeakLags) < maxRequiredIndex
        warning(['%s contains only %d permutation entries, but index %d is ' 'required. Skipping.'], mouseID, numel(allRawPeakLags), maxRequiredIndex);
        continue;
    end

    %% ---------- accepted permutation mask ----------

    if isfield(sess, 'permAccepted') && ...
       ~isempty(sess.permAccepted)

        permAccepted = logical(sess.permAccepted(:)');

        if numel(permAccepted) < maxRequiredIndex
            warning(['%s permAccepted is shorter than the required permutation ' 'range. Validity will be based only on finite peak lags.'], mouseID);
            permAccepted = true(size(allRawPeakLags));
        end
    else
        permAccepted = true(size(allRawPeakLags));
    end

    %% ---------- extract raw independent permutation blocks ----------

    rawPeakLagsH0 = allRawPeakLags(permIndsH0);
    rawPeakLagsH50 = allRawPeakLags(permIndsH50);
    rawPeakLagsHneg50 = allRawPeakLags(permIndsHneg50);

    acceptedH0 = permAccepted(permIndsH0);
    acceptedH50 = permAccepted(permIndsH50);
    acceptedHneg50 = permAccepted(permIndsHneg50);

    %% ---------- remove rejected or invalid permutations ----------

    validMaskH0 = acceptedH0 & isfinite(rawPeakLagsH0);
    validMaskH50 = acceptedH50 & isfinite(rawPeakLagsH50);
    validMaskHneg50 = acceptedHneg50 & isfinite(rawPeakLagsHneg50);
    rawPeakLagsH0Valid = rawPeakLagsH0(validMaskH0);
    rawPeakLagsH50Valid = rawPeakLagsH50(validMaskH50);
    rawPeakLagsHneg50Valid = rawPeakLagsHneg50(validMaskHneg50);

    %% ---------- construct model-specific peak-lag distributions ----------

    peakLagsH0 = rawPeakLagsH0Valid;
    peakLagsH50 = rawPeakLagsH50Valid + lagShiftSec;
    peakLagsHneg50 = rawPeakLagsHneg50Valid - lagShiftSec;
    realPeakLag = sess.peakLag;
    realPeakLags(iSess) = realPeakLag;

    %% ---------- distribution summaries ----------

    [meanH0, medianH0, ci95H0] = summarizeDistribution(peakLagsH0);
    [meanH50, medianH50, ci95H50] = summarizeDistribution(peakLagsH50);
    [meanHneg50, medianHneg50, ci95Hneg50] = summarizeDistribution(peakLagsHneg50);

    %% ---------- empirical compatibility probabilities ----------

    pValH0 = empiricalCenteredTwoSidedProbability(realPeakLag, peakLagsH0);
    pValH50 = empiricalCenteredTwoSidedProbability(realPeakLag, peakLagsH50);
    pValHneg50 = empiricalCenteredTwoSidedProbability(realPeakLag, peakLagsHneg50);

    %% ---------- evidence ratios ----------

    evidenceRatio_H0_over_H50 = safeRatio(pValH0, pValH50);

    evidenceRatio_H0_over_Hneg50 = safeRatio(pValH0, pValHneg50);

    pH0All(iSess) = pValH0;
    pH50All(iSess) = pValH50;
    pHneg50All(iSess) = pValHneg50;

    evidenceH0overH50(iSess) = evidenceRatio_H0_over_H50;
    evidenceH0overHneg50(iSess) = evidenceRatio_H0_over_Hneg50;

    %% ---------- store session output ----------

    R = posthocResults.sessions(iSess);
    R.mouseID = char(mouseID);
   
    if isfield(sess, 'baseSessionName')
        R.baseSessionName = sess.baseSessionName;
    end
   
    if isfield(sess, 'processedDataFolder')
        R.processedDataFolder = sess.processedDataFolder;
    end

    R.realPeakLagSec = realPeakLag;

    if isfield(sess, 'peakCorr')
        R.realPeakCorr = sess.peakCorr;
    end

    R.rawPeakLagsH0 = rawPeakLagsH0;
    R.rawPeakLagsH50 = rawPeakLagsH50;
    R.rawPeakLagsHneg50 = rawPeakLagsHneg50;

    R.peakLagsH0 = peakLagsH0;
    R.peakLagsH50 = peakLagsH50;
    R.peakLagsHneg50 = peakLagsHneg50;

    R.nValidH0 = numel(peakLagsH0);
    R.nValidH50 = numel(peakLagsH50);
    R.nValidHneg50 = numel(peakLagsHneg50);

    R.meanH0 = meanH0;
    R.meanH50 = meanH50;
    R.meanHneg50 = meanHneg50;

    R.medianH0 = medianH0;
    R.medianH50 = medianH50;
    R.medianHneg50 = medianHneg50;

    R.ci95H0 = ci95H0;
    R.ci95H50 = ci95H50;
    R.ci95Hneg50 = ci95Hneg50;

    R.pValH0 = pValH0;
    R.pValH50 = pValH50;
    R.pValHneg50 = pValHneg50;

    R.evidenceRatio_H0_over_H50 = evidenceRatio_H0_over_H50;

    R.evidenceRatio_H0_over_Hneg50 = evidenceRatio_H0_over_Hneg50;

    posthocResults.sessions(iSess) = R;

    %% ---------- print results ----------

    fprintf('real peak lag: %+.3f s\n', realPeakLag);

    fprintf('valid model lags: H0=%d | H+50=%d | H-50=%d\n', numel(peakLagsH0), numel(peakLagsH50), numel(peakLagsHneg50));

    fprintf('model medians: H0=%+.3f | H+50=%+.3f | H-50=%+.3f s\n', medianH0, medianH50, medianHneg50);

    fprintf('probabilities: pH0=%.6f | pH+50=%.6f | pH-50=%.6f\n', pValH0, pValH50, pValHneg50);

    fprintf('evidence ratios: H0/H+50=%.6f | H0/H-50=%.6f\n', evidenceRatio_H0_over_H50, evidenceRatio_H0_over_Hneg50);

    %% ---------- plot distributions ----------

    if makePlots
        plotAnimalPeakLagDistributions(mouseID, realPeakLag, peakLagsH0, peakLagsH50, peakLagsHneg50, ...
            nHistogramBins, pValH0, pValH50, pValHneg50, evidenceRatio_H0_over_H50, evidenceRatio_H0_over_Hneg50);
    end
end

%% ---------------- summary table ----------------

summaryTable = table(animalIDs, realPeakLags, pH0All, pH50All, pHneg50All, evidenceH0overH50, evidenceH0overHneg50, ...
    'VariableNames', {'Animal', 'RealPeakLag_s', 'pValH0', 'pValH50', 'pValHneg50', 'Evidence_H0_over_H50', 'Evidence_H0_over_Hneg50'});

posthocResults.summaryTable = summaryTable;

fprintf('\n================ SUMMARY ================\n');
disp(summaryTable);

%% ---------------- save ----------------

[inputFolder, inputName, ~] = fileparts(savedCrossCorrFile);

outFile = fullfile(inputFolder, inputName + "_PeakLagBayesPosthoc50ms.mat");

save(outFile, 'posthocResults', 'summaryTable', '-v7.3');

fprintf('\nsaved post hoc peak-lag analysis to:\n%s\n', outFile);

end

%% ========================================================================
%% helper: empirical model-centered two-sided probability
%% ========================================================================

function pVal = empiricalCenteredTwoSidedProbability(realValue, modelDistribution)
% Compares the distance of the real value from the model distribution's median with the distances of the model samples from that same median.
% A larger value means the real peak lag is more typical under the model.

modelDistribution = modelDistribution(isfinite(modelDistribution));

if ~isfinite(realValue) || isempty(modelDistribution)
    pVal = NaN;
    return;
end

modelCenter = median(modelDistribution, 'omitnan');

realDeviation = abs(realValue - modelCenter);

nullDeviations = abs(modelDistribution - modelCenter);

pVal = (sum(nullDeviations >= realDeviation) + 1) / (numel(nullDeviations) + 1);

end

%% ========================================================================
%% helper: summarize one model distribution
%% ========================================================================

function [meanValue, medianValue, ci95] = summarizeDistribution(modelDistribution)

modelDistribution = modelDistribution(isfinite(modelDistribution));

if isempty(modelDistribution)
    meanValue = NaN;
    medianValue = NaN;
    ci95 = [NaN NaN];
    return;
end

meanValue = mean(modelDistribution, 'omitnan');
medianValue = median(modelDistribution, 'omitnan');

if numel(modelDistribution) >= 2
    ci95 = prctile(modelDistribution, [2.5 97.5]);
else
    ci95 = [NaN NaN];
end

end

%% ========================================================================
%% helper: safe ratio
%% ========================================================================

function ratio = safeRatio(numerator, denominator)

if ~isfinite(numerator) || ~isfinite(denominator) || denominator == 0

    ratio = NaN;
else
    ratio = numerator / denominator;
end

end

%% ========================================================================
%% helper: plot one animal
%% ========================================================================

function plotAnimalPeakLagDistributions( ...
    mouseID, ...
    realPeakLag, ...
    peakLagsH0, ...
    peakLagsH50, ...
    peakLagsHneg50, ...
    nHistogramBins, ...
    pValH0, ...
    pValH50, ...
    pValHneg50, ...
    evidenceRatioH0H50, ...
    evidenceRatioH0Hneg50)

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
    xPadding = 0.05 * (xMax - xMin);
end

commonXLimits = ...
    [xMin - xPadding, xMax + xPadding];

commonEdges = linspace( ...
    commonXLimits(1), ...
    commonXLimits(2), ...
    nHistogramBins + 1);

figure( ...
    'Name', sprintf('%s Population Peak-Lag Models', mouseID), ...
    'Color', 'w');

tileLayout = tiledlayout( ...
    1, ...
    3, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

title(tileLayout, ...
    sprintf('%s Population Peak-Lag Model Comparison', mouseID), ...
    'FontSize', 17);

%% ---------- H0 ----------

nexttile;
hold on;

histogram( ...
    peakLagsH0, ...
    'BinEdges', commonEdges, ...
    'FaceColor', [0.3 0.6 0.8], ...
    'EdgeColor', 'none');

xline(realPeakLag, 'r-', 'LineWidth', 2);

xlim(commonXLimits);

xlabel('Peak Lag (s)');
ylabel('Count');

title(sprintf('H0 | p = %.3f', pValH0));

box off;
set(gca, 'FontSize', 13, 'TickDir', 'out');

%% ---------- H+50 ----------

nexttile;
hold on;

histogram( ...
    peakLagsH50, ...
    'BinEdges', commonEdges, ...
    'FaceColor', [0.3 0.6 0.8], ...
    'EdgeColor', 'none');

xline(realPeakLag, 'r-', 'LineWidth', 2);

xlim(commonXLimits);

xlabel('Peak Lag (s)');
ylabel('Count');

title(sprintf( ...
    'H+50 | p = %.3f | H0/H+50 = %.3f', ...
    pValH50, ...
    evidenceRatioH0H50));

box off;
set(gca, 'FontSize', 13, 'TickDir', 'out');

%% ---------- H-50 ----------

nexttile;
hold on;

histogram( ...
    peakLagsHneg50, ...
    'BinEdges', commonEdges, ...
    'FaceColor', [0.3 0.6 0.8], ...
    'EdgeColor', 'none');

xline(realPeakLag, 'r-', 'LineWidth', 2);

xlim(commonXLimits);

xlabel('Peak Lag (s)');
ylabel('Count');

title(sprintf( ...
    'H-50 | p = %.3f | H0/H-50 = %.3f', ...
    pValHneg50, ...
    evidenceRatioH0Hneg50));

box off;
set(gca, 'FontSize', 13, 'TickDir', 'out');

end
