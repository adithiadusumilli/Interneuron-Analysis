function computeChunkedPopAvgPeakLagBayesPosthoc50ms(savedCrossCorrFile)
% Post hoc population-averaged peak-lag model comparison, CHUNKED.

% This analysis uses PEAK LAG itself as the model-comparison metric.
% It follows the same statistical logic as the validated no-chunk version.

% Expected input:
%   concatCrossCorr_trialavg_chunked_popavg_ALL_SESSIONS.mat

% The input file must contain:
%   allSessions.sessions(i).mouseID
%   allSessions.sessions(i).real_peakLagSec
%   allSessions.sessions(i).perm_peakLagSec

% Permutation blocks:
%   1:100   = H0 raw peak lags
%   101:200 = H+50 raw peak lags
%   201:300 = H-50 raw peak lags

% Post hoc transformations:
%   H0 peak lags    = raw permutations 1:100
%   H+50 peak lags  = raw permutations 101:200 + 0.050 s
%   H-50 peak lags  = raw permutations 201:300 - 0.050 s

% Real peak lag is compared directly with each model-specific peak-lag
% distribution using one-tailed empirical probabilities with a +1 correction.

% H0 versus H+50 uses the LEFT tail for both distributions:
%   pH0_left = (# H0 peak lags <= real peak lag + 1) / (N_H0 + 1)
%   pH50     = (# H+50 peak lags <= real peak lag + 1) / (N_H50 + 1)

% H0 versus H-50 uses the RIGHT tail for both distributions:
%   pH0_right = (# H0 peak lags >= real peak lag + 1) / (N_H0 + 1)
%   pHneg50   = (# H-50 peak lags >= real peak lag + 1) / (N_Hneg50 + 1)

% Evidence ratios:
%   H0/H+50 = pH0_left  / pH50
%   H0/H-50 = pH0_right / pHneg50

% Interpretation:
%   ratio > 1  favors H0
%   ratio < 1  favors the shifted model
%   ratio ~ 1  indicates little distinction

% Example:
% computeChunkedPopAvgPeakLagBayesPosthoc50ms("C:\Users\mirilab\Documents\GlobusTransfer\concatCrossCorr_trialavg_chunked_popavg_ALL_SESSIONS.mat")

arguments
    savedCrossCorrFile (1,1) string
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
    error('Input file not found:\n%s', savedCrossCorrFile);
end

S = load(savedCrossCorrFile, 'allSessions');

if ~isfield(S, 'allSessions')
    error('Input file does not contain allSessions.');
end

allSessions = S.allSessions;

if ~isfield(allSessions, 'sessions') || isempty(allSessions.sessions)
    error('allSessions.sessions is missing or empty.');
end

sessions = allSessions.sessions;
nSess = numel(sessions);

%% ---------------- initialize output ----------------

posthocResults = struct();
posthocResults.sourceFile = savedCrossCorrFile;
posthocResults.lagShiftSec = lagShiftSec;
posthocResults.metric = 'chunked population peak lag';
posthocResults.noSignFlipApplied = true;

posthocResults.analysisDescription = [ ...
    'The real chunked population-average peak lag was compared with three ' ...
    'independent label-permutation peak-lag distributions. Permutations ' ...
    '1:100 were used unshifted for H0, permutations 101:200 were shifted ' ...
    'by +50 ms for H+50, and permutations 201:300 were shifted by -50 ms ' ...
    'for H-50. No sign flipping was applied.'];

posthocResults.probabilityDefinition = [ ...
    'Empirical one-tailed probabilities compare the real peak lag directly ' ...
    'with each model distribution. H0 versus H+50 uses left-tail ' ...
    'probabilities (permutation lag <= real lag), whereas H0 versus H-50 ' ...
    'uses right-tail probabilities (permutation lag >= real lag). A +1 ' ...
    'correction is applied to numerator and denominator.'];

posthocResults.permIndsH0 = permIndsH0;
posthocResults.permIndsH50 = permIndsH50;
posthocResults.permIndsHneg50 = permIndsHneg50;

emptySessionResult = struct( ...
    'mouseID', '', ...
    'sessionTag', '', ...
    'baseDir', '', ...
    'realPeakLagSec', NaN, ...
    'rawPeakLagsH0', [], ...
    'rawPeakLagsH50', [], ...
    'rawPeakLagsHneg50', [], ...
    'peakLagsH0', [], ...
    'peakLagsH50', [], ...
    'peakLagsHneg50', [], ...
    'nValidH0', 0, ...
    'nValidH50', 0, ...
    'nValidHneg50', 0, ...
    'meanH0', NaN, ...
    'meanH50', NaN, ...
    'meanHneg50', NaN, ...
    'medianH0', NaN, ...
    'medianH50', NaN, ...
    'medianHneg50', NaN, ...
    'ci95H0', [NaN NaN], ...
    'ci95H50', [NaN NaN], ...
    'ci95Hneg50', [NaN NaN], ...
    'pValH0_left', NaN, ...
    'pValH0_right', NaN, ...
    'pValH50', NaN, ...
    'pValHneg50', NaN, ...
    'evidenceRatio_H0_over_H50', NaN, ...
    'evidenceRatio_H0_over_Hneg50', NaN);

posthocResults.sessions = repmat(emptySessionResult, nSess, 1);

%% ---------------- summary arrays ----------------

animalIDs = strings(nSess,1);
sessionTags = strings(nSess,1);
realPeakLags = nan(nSess,1);
nValidH0All = zeros(nSess,1);
nValidH50All = zeros(nSess,1);
nValidHneg50All = zeros(nSess,1);

pH0LeftAll = nan(nSess,1);
pH0RightAll = nan(nSess,1);
pH50All = nan(nSess,1);
pHneg50All = nan(nSess,1);

evidenceH0overH50 = nan(nSess,1);
evidenceH0overHneg50 = nan(nSess,1);

%% ========================================================================
%% loop through sessions
%% ========================================================================

for iSess = 1:nSess

    sess = sessions(iSess);

    %% ---------- identifiers ----------

    if isfield(sess, 'mouseID') && ~isempty(sess.mouseID)
        mouseID = string(sess.mouseID);
    else
        mouseID = "Animal_" + iSess;
    end

    if isfield(sess, 'sessionTag') && ~isempty(sess.sessionTag)
        sessionTag = string(sess.sessionTag);
    else
        sessionTag = mouseID;
    end

    animalIDs(iSess) = mouseID;
    sessionTags(iSess) = sessionTag;

    fprintf('\n========================================\n');
    fprintf('processing %s (%s)\n', mouseID, sessionTag);
    fprintf('========================================\n');

    %% ---------- validate required fields ----------

    if ~isfield(sess, 'real_peakLagSec') || ...
       isempty(sess.real_peakLagSec) || ...
       ~isscalar(sess.real_peakLagSec) || ...
       ~isfinite(sess.real_peakLagSec)

        warning('%s has no valid real_peakLagSec. Skipping.', sessionTag);
        continue;
    end

    if ~isfield(sess, 'perm_peakLagSec') || isempty(sess.perm_peakLagSec)
        warning('%s is missing perm_peakLagSec. Skipping.', sessionTag);
        continue;
    end

    allRawPeakLags = double(sess.perm_peakLagSec(:)');
    maxRequiredIndex = max([permIndsH0, permIndsH50, permIndsHneg50]);

    if numel(allRawPeakLags) < maxRequiredIndex
        warning(['%s contains only %d permutation peak lags, but index %d ' ...
                 'is required. Skipping.'], ...
                 sessionTag, numel(allRawPeakLags), maxRequiredIndex);
        continue;
    end

    %% ---------- extract independent permutation blocks ----------

    rawPeakLagsH0 = allRawPeakLags(permIndsH0);
    rawPeakLagsH50 = allRawPeakLags(permIndsH50);
    rawPeakLagsHneg50 = allRawPeakLags(permIndsHneg50);

    %% ---------- remove invalid entries only ----------

    rawPeakLagsH0Valid = rawPeakLagsH0(isfinite(rawPeakLagsH0));
    rawPeakLagsH50Valid = rawPeakLagsH50(isfinite(rawPeakLagsH50));
    rawPeakLagsHneg50Valid = rawPeakLagsHneg50(isfinite(rawPeakLagsHneg50));

    if isempty(rawPeakLagsH0Valid) || ...
       isempty(rawPeakLagsH50Valid) || ...
       isempty(rawPeakLagsHneg50Valid)

        warning('%s has an empty model distribution after removing invalid values.', sessionTag);
        continue;
    end

    %% ---------- construct model-specific peak-lag distributions ----------

    peakLagsH0 = rawPeakLagsH0Valid;
    peakLagsH50 = rawPeakLagsH50Valid + lagShiftSec;
    peakLagsHneg50 = rawPeakLagsHneg50Valid - lagShiftSec;

    realPeakLag = double(sess.real_peakLagSec);
    realPeakLags(iSess) = realPeakLag;

    %% ---------- distribution summaries ----------

    [meanH0, medianH0, ci95H0] = summarizeDistribution(peakLagsH0);
    [meanH50, medianH50, ci95H50] = summarizeDistribution(peakLagsH50);
    [meanHneg50, medianHneg50, ci95Hneg50] = summarizeDistribution(peakLagsHneg50);

    %% ---------- one-tailed empirical probabilities ----------

    % H0 versus H+50: left tail for both model distributions.
    pValH0_left = empiricalOneTailedProbability( ...
        realPeakLag, peakLagsH0, "left");

    pValH50 = empiricalOneTailedProbability( ...
        realPeakLag, peakLagsH50, "left");

    % H0 versus H-50: right tail for both model distributions.
    pValH0_right = empiricalOneTailedProbability( ...
        realPeakLag, peakLagsH0, "right");

    pValHneg50 = empiricalOneTailedProbability( ...
        realPeakLag, peakLagsHneg50, "right");

    %% ---------- evidence ratios ----------

    evidenceRatio_H0_over_H50 = safeRatio(pValH0_left, pValH50);
    evidenceRatio_H0_over_Hneg50 = safeRatio(pValH0_right, pValHneg50);

    %% ---------- fill summary arrays ----------

    nValidH0All(iSess) = numel(peakLagsH0);
    nValidH50All(iSess) = numel(peakLagsH50);
    nValidHneg50All(iSess) = numel(peakLagsHneg50);

    pH0LeftAll(iSess) = pValH0_left;
    pH0RightAll(iSess) = pValH0_right;
    pH50All(iSess) = pValH50;
    pHneg50All(iSess) = pValHneg50;

    evidenceH0overH50(iSess) = evidenceRatio_H0_over_H50;
    evidenceH0overHneg50(iSess) = evidenceRatio_H0_over_Hneg50;

    %% ---------- store session output ----------

    R = emptySessionResult;
    R.mouseID = char(mouseID);
    R.sessionTag = char(sessionTag);

    if isfield(sess, 'baseDir') && ~isempty(sess.baseDir)
        R.baseDir = sess.baseDir;
    end

    R.realPeakLagSec = realPeakLag;

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

    R.pValH0_left = pValH0_left;
    R.pValH0_right = pValH0_right;
    R.pValH50 = pValH50;
    R.pValHneg50 = pValHneg50;

    R.evidenceRatio_H0_over_H50 = evidenceRatio_H0_over_H50;
    R.evidenceRatio_H0_over_Hneg50 = evidenceRatio_H0_over_Hneg50;

    posthocResults.sessions(iSess) = R;

    %% ---------- print results ----------

    fprintf('real peak lag: %+.3f s\n', realPeakLag);

    fprintf('valid model lags: H0=%d | H+50=%d | H-50=%d\n', ...
        numel(peakLagsH0), numel(peakLagsH50), numel(peakLagsHneg50));

    fprintf('model medians: H0=%+.3f | H+50=%+.3f | H-50=%+.3f s\n', ...
        medianH0, medianH50, medianHneg50);

    fprintf('left-tail probabilities:  pH0=%.6f | pH+50=%.6f\n', ...
        pValH0_left, pValH50);

    fprintf('right-tail probabilities: pH0=%.6f | pH-50=%.6f\n', ...
        pValH0_right, pValHneg50);

    fprintf('evidence ratios: H0/H+50=%.6f | H0/H-50=%.6f\n', ...
        evidenceRatio_H0_over_H50, evidenceRatio_H0_over_Hneg50);

    %% ---------- plot distributions ----------

    if makePlots
        plotAnimalPeakLagDistributions( ...
            mouseID, ...
            sessionTag, ...
            realPeakLag, ...
            peakLagsH0, ...
            peakLagsH50, ...
            peakLagsHneg50, ...
            nHistogramBins, ...
            pValH0_left, ...
            pValH0_right, ...
            pValH50, ...
            pValHneg50, ...
            evidenceRatio_H0_over_H50, ...
            evidenceRatio_H0_over_Hneg50);
    end
end

%% ---------------- summary table ----------------

summaryTable = table( ...
    animalIDs, ...
    sessionTags, ...
    realPeakLags, ...
    nValidH0All, ...
    nValidH50All, ...
    nValidHneg50All, ...
    pH0LeftAll, ...
    pH0RightAll, ...
    pH50All, ...
    pHneg50All, ...
    evidenceH0overH50, ...
    evidenceH0overHneg50, ...
    'VariableNames', { ...
        'Animal', ...
        'SessionTag', ...
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

fprintf('\n================ SUMMARY ================\n');
disp(summaryTable);

%% ---------------- save ----------------

[inputFolder, ~, ~] = fileparts(savedCrossCorrFile);

outFile = fullfile( ...
    inputFolder, ...
    "concatCrossCorr_trialavg_chunked_popavg_ALL_SESSIONS_PeakLagBayesPosthoc50ms.mat");

save(outFile, 'posthocResults', 'summaryTable', '-v7.3');

fprintf('\nsaved post hoc chunked population peak-lag analysis to:\n%s\n', ...
    outFile);

end

%% ========================================================================
%% helper: empirical one-tailed probability
%% ========================================================================

function pVal = empiricalOneTailedProbability(realValue, modelDistribution, tail)
% Compare the real peak lag directly with one model distribution.
%
% left:  P(X <= realValue)
% right: P(X >= realValue)
%
% The +1 correction prevents a zero empirical probability.

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

pVal = (nExtreme + 1) / (numel(modelDistribution) + 1);

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

if ~isfinite(numerator) || ...
   ~isfinite(denominator) || ...
   denominator == 0

    ratio = NaN;
else
    ratio = numerator / denominator;
end

end

%% ========================================================================
%% helper: plot one session
%% ========================================================================

function plotAnimalPeakLagDistributions( ...
    mouseID, ...
    sessionTag, ...
    realPeakLag, ...
    peakLagsH0, ...
    peakLagsH50, ...
    peakLagsHneg50, ...
    nHistogramBins, ...
    pValH0Left, ...
    pValH0Right, ...
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

commonXLimits = [xMin - xPadding, xMax + xPadding];

commonEdges = linspace( ...
    commonXLimits(1), ...
    commonXLimits(2), ...
    nHistogramBins + 1);

figure( ...
    'Name', sprintf('%s Chunked Population Peak-Lag Models', mouseID), ...
    'Color', 'w');

tileLayout = tiledlayout( ...
    1, ...
    3, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

if string(sessionTag) == string(mouseID)
    overallTitle = sprintf( ...
        '%s Chunked Population Peak-Lag Model Comparison', mouseID);
else
    overallTitle = sprintf( ...
        '%s | %s | Chunked Population Peak-Lag Model Comparison', ...
        mouseID, sessionTag);
end

title(tileLayout, overallTitle, 'FontSize', 17);

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

title(sprintf( ...
    'H0 | p_{left} = %.3f | p_{right} = %.3f', ...
    pValH0Left, pValH0Right));

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
    'H+50 | p_{left} = %.3f | H0/H+50 = %.3f', ...
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
    'H-50 | p_{right} = %.3f | H0/H-50 = %.3f', ...
    pValHneg50, ...
    evidenceRatioH0Hneg50));

box off;
set(gca, 'FontSize', 13, 'TickDir', 'out');

end
