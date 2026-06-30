function classifyCortexAndStriatum_getMouseDataNames(mouseIDs, baseSessionNames, probeRegions, outFile)
% classifies cortex and striatum neurons using waveform width + GMM
% this version uses David's getMouseDataNames cuz new animals have some diff file naming

% Classification method follows Miri et al. (2017, Neuron):
%   1. Pool all neurons across all animals into one distribution per region
%   2. Fit a 2-component GMM to the pooled widths (0 to 0.8 ms range)
%   3. Find two asymmetric width boundaries using CDF-based misclassification
%      rate equations, at a prescribed error rate (default 1%)
%      - lowerBound: below this = narrow (interneuron / FSI)
%      - upperBound: above this = wide (pyramidal / MSN)
%      - between bounds = unclassified (too close to overlap region)
%   4. Apply those global boundaries to classify neurons per animal

% Outputs saved to AA_classifications.mat include:
%   - classifications: {nFiles x 2} cell of per-neuron labels (0=wide, 1=narrow, NaN=wrong region OR unclassified)
%   - classificationsWithGap: same but unclassified gap neurons marked -1 instead of NaN (use to identify/count them)
%   - fileWidths: {nFiles x 2} cell of per-neuron widths
%   - intersectionPointCortex / intersectionPointStriat: GMM intersection (for reference)
%   - lowerBoundCortex / upperBoundCortex: misclassification boundaries
%   - lowerBoundStriat / upperBoundStriat
%   - gmmParamsCortex / gmmParamsStriat: struct with mu, sigma, weights

% order for new AA_classifications:
%   1 D026
%   2 D020
%   3 D024
%   4 D043
%   5 D050
%   6 D054

if nargin < 4 || isempty(outFile)
    outFile = 'AA_classifications.mat';
end

% prescribed misclassification rate (fraction)
misclassThresh = 0.01;  % 1%

%% make sure MATLAB can find getMouseDataNames.m
addpath('C:\Users\mirilab\Documents\GlobusTransfer');

consolidatedDataFolder = 'X:\David\AnalysesData';

numFiles = numel(mouseIDs);

if numel(baseSessionNames) ~= numFiles || numel(probeRegions) ~= numFiles
    error('mouseIDs, baseSessionNames, and probeRegions must all have the same length');
end

classifications        = cell(numFiles,2); % col 1 = cortex, col 2 = striatum
classificationsWithGap = cell(numFiles,2); % same but gap neurons marked -1 instead of NaN

fileWidths = cell(numFiles,2); % col 1 = cortex widths, col 2 = striatum widths

neuronDataStructFiles = cell(numFiles,1);
firingRatesFiles      = cell(numFiles,1);

%% ---------------- get all file paths from David's helper ----------------
for fileIndex = 1:numFiles

    dataNames = getMouseDataNames( ...
        mouseIDs{fileIndex}, ...
        baseSessionNames{fileIndex}, ...
        probeRegions{fileIndex});

    neuronDataStructFiles{fileIndex} = dataNames.neuronDataStruct;
    firingRatesFiles{fileIndex}      = dataNames.NeuralFiringRates1msBins10msGauss;

    fprintf('\nfile %d: %s\n', fileIndex, mouseIDs{fileIndex});
    fprintf('neuronDataStruct: %s\n', neuronDataStructFiles{fileIndex});
    fprintf('firing rates:      %s\n', firingRatesFiles{fileIndex});
end

%% ---------------- PART 1: calculate spike widths ----------------
% width = time between waveform trough and subsequent peak on the biggest channel
% (trough-to-peak, following PI's paper convention)

for fileIndex = 1:numFiles

    load(neuronDataStructFiles{fileIndex}, 'neuronDataStruct');
    load(firingRatesFiles{fileIndex}, 'cortexInds', 'striatumInds');

    %% cortex widths
    cortexWidths = nan(1, numel(cortexInds));

    for i = 1:numel(cortexInds)
        waveform   = neuronDataStruct(cortexInds(i)).waveforms;
        biggestChan = neuronDataStruct(cortexInds(i)).biggestChan;

        ap = waveform(:, biggestChan);

        [~, mn] = min(ap); % trough
        [~, mx] = max(ap(mn:end)); % subsequent peak after trough
        mx = mx + mn - 1; % adjust index to full waveform

        cortexWidths(i) = (mx - mn) / 30000;  % seconds
    end

    fileWidths{fileIndex, 1} = cortexWidths;

    %% striatum widths
    striatWidths = nan(1, numel(striatumInds));

    for i = 1:numel(striatumInds)
        waveform    = neuronDataStruct(striatumInds(i)).waveforms;
        biggestChan = neuronDataStruct(striatumInds(i)).biggestChan;

        ap = waveform(:, biggestChan);

        [~, mn] = min(ap);
        [~, mx] = max(ap(mn:end));
        mx = mx + mn - 1;

        striatWidths(i) = (mx - mn) / 30000;  % seconds
    end

    fileWidths{fileIndex, 2} = striatWidths;
end

%% ---------------- PART 2: pool widths and fit a single global GMM ----------------
% Per PI's method: combine neurons from ALL animals, fit one GMM,
% then derive global boundaries to apply uniformly across animals.
% Only use widths <= 0.8 ms for the fit (as in the paper).

maxWidth = 0.8e-3;  % 0.8 ms in seconds
numComponents = 2;

%% cortex: pool and fit
allCortexWidths = [];
for fileIndex = 1:numFiles
    w = fileWidths{fileIndex, 1};
    w = w(~isnan(w) & isfinite(w) & w <= maxWidth);
    allCortexWidths = [allCortexWidths, w]; %#ok<AGROW>
end

fprintf('\nFitting cortex GMM on %d pooled neurons (across all animals)...\n', numel(allCortexWidths));
gmCortex = fitgmdist(allCortexWidths', numComponents, 'Options', statset('MaxIter', 500));

[muC, orderC]  = sort(gmCortex.mu, 'ascend');
sigC           = squeeze(sqrt(gmCortex.Sigma));
sigC           = sigC(orderC);
wtC            = gmCortex.ComponentProportion(orderC);

gmmParamsCortex = struct('mu', muC, 'sigma', sigC, 'weights', wtC);

intersectionPointCortex = calculateIntersectionPoint(muC, sigC);

[lowerBoundCortex, upperBoundCortex] = calculateMisclassificationBounds( ...
    muC, sigC, misclassThresh);

fprintf('Cortex GMM:  mu=[%.4f  %.4f] ms, sigma=[%.4f  %.4f] ms\n', ...
    muC(1)*1e3, muC(2)*1e3, sigC(1)*1e3, sigC(2)*1e3);
fprintf('Cortex bounds: lower=%.4f ms, upper=%.4f ms\n', ...
    lowerBoundCortex*1e3, upperBoundCortex*1e3);

%% striatum: pool and fit
allStriatWidths = [];
for fileIndex = 1:numFiles
    w = fileWidths{fileIndex, 2};
    w = w(~isnan(w) & isfinite(w) & w <= maxWidth);
    allStriatWidths = [allStriatWidths, w]; %#ok<AGROW>
end

fprintf('\nFitting striatum GMM on %d pooled neurons (across all animals)...\n', numel(allStriatWidths));
gmStriat = fitgmdist(allStriatWidths', numComponents, 'Options', statset('MaxIter', 500));

[muS, orderS]  = sort(gmStriat.mu, 'ascend');
sigS           = squeeze(sqrt(gmStriat.Sigma));
sigS           = sigS(orderS);
wtS            = gmStriat.ComponentProportion(orderS);

gmmParamsStriat = struct('mu', muS, 'sigma', sigS, 'weights', wtS);

intersectionPointStriat = calculateIntersectionPoint(muS, sigS);

[lowerBoundStriat, upperBoundStriat] = calculateMisclassificationBounds( ...
    muS, sigS, misclassThresh);

fprintf('Striatum GMM: mu=[%.4f  %.4f] ms, sigma=[%.4f  %.4f] ms\n', ...
    muS(1)*1e3, muS(2)*1e3, sigS(1)*1e3, sigS(2)*1e3);
fprintf('Striatum bounds: lower=%.4f ms, upper=%.4f ms\n', ...
    lowerBoundStriat*1e3, upperBoundStriat*1e3);

%% ---------------- PART 3: plot pooled distributions with global bounds ----------------

%% cortex pooled plot
figure;
h = histogram(allCortexWidths * 1e3, 'BinWidth', (1/20000)*1e3, ...
    'EdgeColor', 'black', 'FaceColor', 'blue');
hold on;

x = linspace(min(allCortexWidths), max(allCortexWidths), 1000);
xMs = x * 1e3;

yNar = pdf('Normal', x, muC(1), sigC(1)) * numel(allCortexWidths) * h.BinWidth * 1e-3 * wtC(1);
yWid = pdf('Normal', x, muC(2), sigC(2)) * numel(allCortexWidths) * h.BinWidth * 1e-3 * wtC(2);

plot(xMs, yWid, 'r', 'LineWidth', 2);
plot(xMs, yNar, 'g', 'LineWidth', 2);
xline(lowerBoundCortex * 1e3, 'k--', 'LineWidth', 1.5);
xline(upperBoundCortex * 1e3, 'k:', 'LineWidth', 1.5);

title('ALL ANIMALS Cortex Spike Widths — Global GMM Boundaries');
legend({'Widths','Pyramidal (wide)','Interneuron (narrow)','Lower bound','Upper bound'});
xlabel('Trough-to-peak width (ms)');
ylabel('Count');
box off;

%% striatum pooled plot
figure;
h = histogram(allStriatWidths * 1e3, 'BinWidth', (1/20000)*1e3, ...
    'EdgeColor', 'black', 'FaceColor', 'magenta');
hold on;

x = linspace(min(allStriatWidths), max(allStriatWidths), 1000);
xMs = x * 1e3;

yNar = pdf('Normal', x, muS(1), sigS(1)) * numel(allStriatWidths) * h.BinWidth * 1e-3 * wtS(1);
yWid = pdf('Normal', x, muS(2), sigS(2)) * numel(allStriatWidths) * h.BinWidth * 1e-3 * wtS(2);

plot(xMs, yWid, 'r', 'LineWidth', 2);
plot(xMs, yNar, 'g', 'LineWidth', 2);
xline(lowerBoundStriat * 1e3, 'k--', 'LineWidth', 1.5);
xline(upperBoundStriat * 1e3, 'k:', 'LineWidth', 1.5);

title('ALL ANIMALS Striatum Spike Widths — Global GMM Boundaries');
legend({'Widths','Wide (MSN)','Narrow (FSI)','Lower bound','Upper bound'});
xlabel('Trough-to-peak width (ms)');
ylabel('Count');
box off;

%% per-animal plots (for inspection) using global bounds
for fileIndex = 1:numFiles
    animalID = mouseIDs{fileIndex};

    %% per-animal cortex
    cw = fileWidths{fileIndex, 1};
    cw = cw(~isnan(cw) & isfinite(cw));

    figure;
    h = histogram(cw * 1e3, 'BinWidth', (1/20000)*1e3, ...
        'EdgeColor', 'black', 'FaceColor', 'blue');
    hold on;
    xline(lowerBoundCortex * 1e3, 'k--', 'LineWidth', 1.5);
    xline(upperBoundCortex * 1e3, 'k:', 'LineWidth', 1.5);
    title(sprintf('%s Cortex — Global Boundaries Applied', animalID));
    legend({'Widths','Lower bound (narrow cutoff)','Upper bound (wide cutoff)'});
    xlabel('Trough-to-peak width (ms)');
    ylabel('Count');
    box off;

    %% per-animal striatum
    sw = fileWidths{fileIndex, 2};
    sw = sw(~isnan(sw) & isfinite(sw));

    figure;
    h = histogram(sw * 1e3, 'BinWidth', (1/20000)*1e3, ...
        'EdgeColor', 'black', 'FaceColor', 'magenta');
    hold on;
    xline(lowerBoundStriat * 1e3, 'k--', 'LineWidth', 1.5);
    xline(upperBoundStriat * 1e3, 'k:', 'LineWidth', 1.5);
    title(sprintf('%s Striatum — Global Boundaries Applied', animalID));
    legend({'Widths','Lower bound (narrow cutoff)','Upper bound (wide cutoff)'});
    xlabel('Trough-to-peak width (ms)');
    ylabel('Count');
    box off;
end

%% ---------------- PART 4: assign labels using global boundaries ----------------
% label convention for `classifications` (drop-in replacement for downstream code):
%   0   = wide waveform  (pyramidal / MSN)
%   1   = narrow waveform (interneuron / FSI)
%   NaN = not in that region  OR  unclassified (in gap between bounds)
%
% This means all existing downstream code that uses == 0, == 1, or ~isnan()
% will automatically exclude gap neurons without any changes.
%
% `classificationsWithGap` is identical except gap neurons are marked -1
% instead of NaN, so you can identify and count them if needed.

for fileIndex = 1:numFiles

    load(neuronDataStructFiles{fileIndex}, 'neuronDataStruct');
    load(firingRatesFiles{fileIndex}, 'cortexInds', 'striatumInds');

    cortexLabels     = nan(1, numel(neuronDataStruct));
    striatLabels     = nan(1, numel(neuronDataStruct));

    cortexLabelsGap  = nan(1, numel(neuronDataStruct));
    striatLabelsGap  = nan(1, numel(neuronDataStruct));

    %% cortex labels
    for i = 1:numel(cortexInds)
        w = fileWidths{fileIndex, 1}(i);

        if w < lowerBoundCortex
            label = 1;   % narrow = interneuron
        elseif w > upperBoundCortex
            label = 0;   % wide = pyramidal
        else
            label = NaN; % unclassified (in gap between bounds)
        end

        cortexLabels(cortexInds(i))    = label;
        cortexLabelsGap(cortexInds(i)) = label;
        if isnan(label)
            cortexLabelsGap(cortexInds(i)) = -1;  % mark gap neurons explicitly
        end
    end

    %% striatum labels
    for i = 1:numel(striatumInds)
        w = fileWidths{fileIndex, 2}(i);

        if w < lowerBoundStriat
            label = 1;   % narrow = FSI
        elseif w > upperBoundStriat
            label = 0;   % wide = MSN
        else
            label = NaN; % unclassified
        end

        striatLabels(striatumInds(i))    = label;
        striatLabelsGap(striatumInds(i)) = label;
        if isnan(label)
            striatLabelsGap(striatumInds(i)) = -1;
        end
    end

    classifications{fileIndex, 1}        = cortexLabels;
    classifications{fileIndex, 2}        = striatLabels;
    classificationsWithGap{fileIndex, 1} = cortexLabelsGap;
    classificationsWithGap{fileIndex, 2} = striatLabelsGap;
end

%% ---------------- PART 5: save everything ----------------
outPath = fullfile(consolidatedDataFolder, outFile);

% human-readable boundary summary (so you never have to rerun to check cutoffs)
boundaryInfo.lowerBoundCortex_ms  = lowerBoundCortex * 1000;
boundaryInfo.upperBoundCortex_ms  = upperBoundCortex * 1000;
boundaryInfo.lowerBoundStriat_ms  = lowerBoundStriat * 1000;
boundaryInfo.upperBoundStriat_ms  = upperBoundStriat * 1000;
boundaryInfo.misclassThresh       = misclassThresh;
boundaryInfo.note = 'lower: below this = narrow; upper: above this = wide; between = unclassified';

save(outPath, ...
    'classifications', ...
    'classificationsWithGap', ...
    'mouseIDs', ...
    'baseSessionNames', ...
    'probeRegions', ...
    'neuronDataStructFiles', ...
    'firingRatesFiles', ...
    'fileWidths', ...
    'intersectionPointCortex', ...
    'intersectionPointStriat', ...
    'lowerBoundCortex', ...
    'upperBoundCortex', ...
    'lowerBoundStriat', ...
    'upperBoundStriat', ...
    'boundaryInfo', ...
    'gmmParamsCortex', ...
    'gmmParamsStriat', ...
    'misclassThresh', ...
    '-v7.3');

fprintf('\nsaved combined cortex & striatum classifications to:\n%s\n', outPath);
fprintf('(misclassification threshold used: %.1f%%)\n', misclassThresh * 100);

fprintf('\n%-8s  %-6s  %-6s  %-6s  ||  %-6s  %-6s  %-6s\n', ...
    'Animal', 'CxPyr', 'CxInt', 'CxUnc', 'StWid', 'StNar', 'StUnc');
for f = 1:numFiles
    fprintf('%-8s  %-6d  %-6d  %-6d  ||  %-6d  %-6d  %-6d\n', ...
        mouseIDs{f}, ...
        sum(classifications{f,1}==0,   'omitnan'), ...
        sum(classifications{f,1}==1,   'omitnan'), ...
        sum(classificationsWithGap{f,1}==-1), ...
        sum(classifications{f,2}==0,   'omitnan'), ...
        sum(classifications{f,2}==1,   'omitnan'), ...
        sum(classificationsWithGap{f,2}==-1));
end

end

%% ========================================================================
%% helpers
%% ========================================================================

function intersectionPoint = calculateIntersectionPoint(means, stdDevs)
% finds where the two unweighted Gaussians cross (used for reference / plotting)
% matches the PI's paper which does not apply mixture weights to the intersection

gaussPDF = @(x, mu, sig) ...
    (1 / (sig * sqrt(2*pi))) * exp(-(x - mu).^2 / (2 * sig^2));

intersectionEquation = @(x) ...
    gaussPDF(x, means(1), stdDevs(1)) - gaussPDF(x, means(2), stdDevs(2));

intersectionPoint = fzero(intersectionEquation, mean(means));

end


function [lowerBound, upperBound] = calculateMisclassificationBounds(means, stdDevs, thresh)
% Finds two asymmetric classification boundaries matching Miri et al. 2017 (Neuron) exactly.
% Equations from the paper (no mixture weights — assumes equal priors):
%
%   lowerBound b_L : mN(b_L) <= thresh
%     fraction of "narrow"-classified neurons that are actually wide
%     mN(b) = cdfW(b) / (cdfW(b) + cdfN(b))
%
%   upperBound b_U : mW(b_U) <= thresh
%     fraction of "wide"-classified neurons that are actually narrow
%     mW(b) = (1-cdfN(b)) / ((1-cdfN(b)) + (1-cdfW(b)))
%
% means(1), stdDevs(1) = narrow distribution (smaller mean)
% means(2), stdDevs(2) = wide distribution (larger mean)

muN  = means(1);    sigN = stdDevs(1);
muW  = means(2);    sigW = stdDevs(2);

cdfN = @(x) 0.5 * (1 + erf((x - muN) / (sigN * sqrt(2))));
cdfW = @(x) 0.5 * (1 + erf((x - muW) / (sigW * sqrt(2))));

% mN(b): fraction of narrow-classified that are actually wide
mN_eq = @(b) cdfW(b) / (cdfW(b) + cdfN(b)) - thresh;

% mW(b): fraction of wide-classified that are actually narrow
mW_eq = @(b) (1 - cdfN(b)) / ((1 - cdfN(b)) + (1 - cdfW(b))) - thresh;

x0 = mean(means);

bNarrowEq = fzero(mN_eq, x0);  % root of the narrow-side equation
bWideEq   = fzero(mW_eq, x0);  % root of the wide-side equation

% fzero can converge to either side depending on equation shape, so don't
% assume which root is smaller — sort them explicitly to get the true
% lower (narrow cutoff) and upper (wide cutoff) bounds, guaranteeing a
% gap between them whenever one exists
lowerBound = min(bNarrowEq, bWideEq);
upperBound = max(bNarrowEq, bWideEq);

end
