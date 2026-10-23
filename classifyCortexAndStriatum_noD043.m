function classifyCortexAndStriatum_noD043(mouseIDs, baseSessionNames, probeRegions)
% Refit the original pooled cortex/striatum GMMs without D043.
% Keeps AA_classifications.mat untouched; writes AA_classifications(nod043).mat.
% Also saves per-animal comparison tables and changed-neuron details.
% Pass the same six-session input lists as the original function.

% classifies cortex and striatum neurons using waveform width + GMM
% this version uses David's getMouseDataNames cuz new animals have some diff file naming

% Classification method follows Andrew's 2017 Neuron paper:
%   1. Pool all neurons across all animals into 1 dist per region
%   2. Fit a 2-component GMM to the pooled widths (0 to 0.8 ms range)
%   3. Find 2 asymmetric width boundaries using CDF-based misclassification
%      rate equations, at a prescribed error rate (default 1%)
%      - lowerBound: below this = narrow (interneuron / FSI)
%      - upperBound: above this = wide (pyramidal / MSN)
%      - between bounds = unclassified (too close to overlap region)
%   4. Apply those global boundaries to classify neurons per animal

% Outputs im saving to AA_classifications.mat (overwrite) include:
%   - classifications: {nFiles x 2} cell of per-neuron labels (0=wide, 1=narrow, NaN=wrong region OR unclassified)
%   - classificationsWithGap: same but unclassified gap neurons marked -1 instead of NaN (use to identify/count them)
%   - fileWidths: {nFiles x 2} cell of per-neuron widths
%   - intersectionPointCortex / intersectionPointStriat: GMM intersection (for reference)
%   - lowerBoundCortex / upperBoundCortex: misclassification boundaries
%   - lowerBoundStriat / upperBoundStriat
%   - gmmParamsCortex / gmmParamsStriat: struct with mu, sigma, weights

% Output rows: D026, D020, D024, D050, D054 (D043 excluded).
% Comparison is matched by mouse ID and original global neuron index.

if nargin < 3
    error('Supply mouseIDs, baseSessionNames, and probeRegions.');
end
mouseIDs = cellstr(string(mouseIDs(:)));
baseSessionNames = cellstr(string(baseSessionNames(:)));
probeRegions = cellstr(string(probeRegions(:)));

keep = ~strcmpi(mouseIDs,'D043');
if ~any(~keep)
    warning('D043 was not in the input list; using all supplied sessions.');
end
mouseIDs = mouseIDs(keep);
baseSessionNames = baseSessionNames(keep);
probeRegions = probeRegions(keep);

outFile = 'AA_classifications(nod043).mat';
oldPath = fullfile('X:\David\AnalysesData','AA_classifications.mat');
if ~isfile(oldPath)
    error('Original classifications missing: %s',oldPath);
end
old = load(oldPath,'classifications','classificationsWithGap','mouseIDs', ...
    'lowerBoundCortex','upperBoundCortex','lowerBoundStriat','upperBoundStriat');
if ~isfield(old,'mouseIDs')
    error('Old file lacks mouseIDs, needed to match animals correctly.');
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
    firingRatesFiles{fileIndex} = dataNames.NeuralFiringRates1msBins10msGauss;

    fprintf('\nfile %d: %s\n', fileIndex, mouseIDs{fileIndex});
    fprintf('neuronDataStruct: %s\n', neuronDataStructFiles{fileIndex});
    fprintf('firing rates: %s\n', firingRatesFiles{fileIndex});
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
        waveform = neuronDataStruct(cortexInds(i)).waveforms;
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
        waveform = neuronDataStruct(striatumInds(i)).waveforms;
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
sigC = squeeze(sqrt(gmCortex.Sigma));
sigC = sigC(orderC);
wtC = gmCortex.ComponentProportion(orderC);

gmmParamsCortex = struct('mu', muC, 'sigma', sigC, 'weights', wtC);

intersectionPointCortex = calculateIntersectionPoint(muC, sigC);

[lowerBoundCortex, upperBoundCortex] = calculateMisclassificationBounds(muC, sigC, misclassThresh);

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

[lowerBoundStriat, upperBoundStriat] = calculateMisclassificationBounds(muS, sigS, misclassThresh);

fprintf('Striatum GMM: mu=[%.4f  %.4f] ms, sigma=[%.4f  %.4f] ms\n', muS(1)*1e3, muS(2)*1e3, sigS(1)*1e3, sigS(2)*1e3);
fprintf('Striatum bounds: lower=%.4f ms, upper=%.4f ms\n', lowerBoundStriat*1e3, upperBoundStriat*1e3);

%% ---------------- PART 3: plot pooled distributions with global bounds ----------------

%% cortex pooled plot
figure;
h = histogram(allCortexWidths * 1e3, 'BinWidth', (1/20000)*1e3, 'EdgeColor', 'black', 'FaceColor', 'blue');
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

% This means all existing downstream code that uses == 0, == 1, or ~isnan() will automatically exclude gap neurons without any changes.

% `classificationsWithGap` is identical except gap neurons are marked -1
% instead of NaN, so you can identify and count them if needed.

for fileIndex = 1:numFiles

    load(neuronDataStructFiles{fileIndex}, 'neuronDataStruct');
    load(firingRatesFiles{fileIndex}, 'cortexInds', 'striatumInds');

    cortexLabels = nan(1, numel(neuronDataStruct));
    striatLabels = nan(1, numel(neuronDataStruct));

    cortexLabelsGap = nan(1, numel(neuronDataStruct));
    striatLabelsGap = nan(1, numel(neuronDataStruct));

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

        cortexLabels(cortexInds(i)) = label;
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

    classifications{fileIndex, 1} = cortexLabels;
    classifications{fileIndex, 2} = striatLabels;
    classificationsWithGap{fileIndex, 1} = cortexLabelsGap;
    classificationsWithGap{fileIndex, 2} = striatLabelsGap;
end

%% ---------------- PART 5: save everything ----------------
outPath = fullfile(consolidatedDataFolder, outFile);

% human-readable boundary summary (so you never have to rerun to check cutoffs)
boundaryInfo.lowerBoundCortex_ms = lowerBoundCortex * 1000;
boundaryInfo.upperBoundCortex_ms = upperBoundCortex * 1000;
boundaryInfo.lowerBoundStriat_ms = lowerBoundStriat * 1000;
boundaryInfo.upperBoundStriat_ms = upperBoundStriat * 1000;
boundaryInfo.misclassThresh = misclassThresh;
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
        sum(classifications{f,1}==0, 'omitnan'), ...
        sum(classifications{f,1}==1, 'omitnan'), ...
        sum(classificationsWithGap{f,1}==-1), ...
        sum(classifications{f,2}==0, 'omitnan'), ...
        sum(classifications{f,2}==1, 'omitnan'), ...
        sum(classificationsWithGap{f,2}==-1));
end


%% ---------------- PART 6: compare old vs new, by animal and region ----------------
% Compare only units belonging to each region in the NEW firing-rate file.
% This avoids treating out-of-region NaNs as unclassified neurons.
animalCol = strings(0,1);
regionCol = strings(0,1);
nUnitsCol = zeros(0,1);
oldWideCol = zeros(0,1); oldNarrowCol = zeros(0,1); oldGapCol = zeros(0,1);
newWideCol = zeros(0,1); newNarrowCol = zeros(0,1); newGapCol = zeros(0,1);
changedCol = zeros(0,1);
oldClassCol = strings(0,1); newClassCol = strings(0,1);
detailAnimal = strings(0,1); detailRegion = strings(0,1);
detailIndex = zeros(0,1); detailWidthMs = zeros(0,1);
detailOld = strings(0,1); detailNew = strings(0,1);

oldIDs = string(old.mouseIDs(:));
for f = 1:numFiles
    oldRow = find(strcmpi(oldIDs,string(mouseIDs{f})),1);
    if isempty(oldRow)
        error('Mouse %s is absent from original AA_classifications.mat.',mouseIDs{f});
    end
    fr = load(firingRatesFiles{f},'cortexInds','striatumInds');
    for r = 1:2
        if r == 1
            regionName = "Cortex";
            inds = fr.cortexInds(:);
        else
            regionName = "Striatum";
            inds = fr.striatumInds(:);
        end
        oldLabels = old.classifications{oldRow,r}(:);
        newLabels = classifications{f,r}(:);
        if any(inds > numel(oldLabels)) || any(inds > numel(newLabels))
            error('Neuron indexing mismatch for %s %s.',mouseIDs{f},regionName);
        end
        a = oldLabels(inds);
        b = newLabels(inds);
        % 0 = wide, 1 = narrow, NaN = unclassified.
        % If original gap labels are available, use them to distinguish gaps.
        if isfield(old,'classificationsWithGap') && ...
                ~isempty(old.classificationsWithGap{oldRow,r})
            oldGapVec = old.classificationsWithGap{oldRow,r}(:);
            oldIsGap = oldGapVec(inds) == -1;
        else
            oldIsGap = isnan(a);
        end
        newGapVec = classificationsWithGap{f,r}(:);
        newIsGap = newGapVec(inds) == -1;

        oldCode = a; oldCode(oldIsGap) = -1;
        newCode = b; newCode(newIsGap) = -1;
        changed = ~(oldCode == newCode | (isnan(oldCode) & isnan(newCode)));

        animalCol(end+1,1) = string(mouseIDs{f});
        regionCol(end+1,1) = regionName;
        nUnitsCol(end+1,1) = numel(inds);
        oldWideCol(end+1,1) = nnz(a==0);
        oldNarrowCol(end+1,1) = nnz(a==1);
        oldGapCol(end+1,1) = nnz(oldIsGap);
        newWideCol(end+1,1) = nnz(b==0);
        newNarrowCol(end+1,1) = nnz(b==1);
        newGapCol(end+1,1) = nnz(newIsGap);
        changedCol(end+1,1) = nnz(changed);

        oldClassCol(end+1,1) = sprintf('wide %d / narrow %d / gap %d', ...
            nnz(a==0),nnz(a==1),nnz(oldIsGap));
        newClassCol(end+1,1) = sprintf('wide %d / narrow %d / gap %d', ...
            nnz(b==0),nnz(b==1),nnz(newIsGap));

        widths = fileWidths{f,r}(:);
        for j = find(changed(:))'
            detailAnimal(end+1,1) = string(mouseIDs{f});
            detailRegion(end+1,1) = regionName;
            detailIndex(end+1,1) = inds(j);
            detailWidthMs(end+1,1) = widths(j)*1000;
            detailOld(end+1,1) = localClassName(oldCode(j));
            detailNew(end+1,1) = localClassName(newCode(j));
        end
    end
end

comparisonTable = table(animalCol,regionCol,nUnitsCol, ...
    oldWideCol,oldNarrowCol,oldGapCol, ...
    newWideCol,newNarrowCol,newGapCol,changedCol, ...
    'VariableNames',{'MouseID','Region','NUnits', ...
    'OldWide','OldNarrow','OldGap','NewWide','NewNarrow','NewGap','NChanged'});
changedNeuronTable = table(detailAnimal,detailRegion,detailIndex,detailWidthMs, ...
    detailOld,detailNew, ...
    'VariableNames',{'MouseID','Region','GlobalNeuronIndex','Width_ms', ...
    'OldClass','NewClass'});

boundaryComparison = table( ...
    ["Cortex";"Cortex";"Striatum";"Striatum"], ...
    ["Lower";"Upper";"Lower";"Upper"], ...
    1000*[old.lowerBoundCortex;old.upperBoundCortex; ...
          old.lowerBoundStriat;old.upperBoundStriat], ...
    1000*[lowerBoundCortex;upperBoundCortex;lowerBoundStriat;upperBoundStriat], ...
    'VariableNames',{'Region','Boundary','Old_ms','New_ms'});
boundaryComparison.Delta_ms = boundaryComparison.New_ms - boundaryComparison.Old_ms;

fprintf('\n=== OLD vs NEW (without D043) ===\n');
disp(comparisonTable);
fprintf('\n=== GMM boundary changes ===\n');
disp(boundaryComparison);
fprintf('\nTotal changed neuron labels: %d\n',height(changedNeuronTable));
disp(changedNeuronTable);

writetable(comparisonTable,fullfile(consolidatedDataFolder, ...
    'AA_classifications_comparison(nod043).csv'));
writetable(changedNeuronTable,fullfile(consolidatedDataFolder, ...
    'AA_classifications_changedNeurons(nod043).csv'));
writetable(boundaryComparison,fullfile(consolidatedDataFolder, ...
    'AA_classifications_boundaries(nod043).csv'));
save(fullfile(consolidatedDataFolder,'AA_classifications_comparison(nod043).mat'), ...
    'comparisonTable','changedNeuronTable','boundaryComparison','-v7.3');

fprintf('\nOriginal preserved: %s\n',oldPath);
fprintf('New classifications: %s\n',outPath);
fprintf('Comparison tables saved in: %s\n',consolidatedDataFolder);

end

function name = localClassName(x)
if x == 0
    name = "Wide";
elseif x == 1
    name = "Narrow";
elseif x == -1
    name = "Gap";
else
    name = "NaN";
end
end

%% ========================================================================
%% helpers
%% ========================================================================

function intersectionPoint = calculateIntersectionPoint(means, stdDevs)
% finds where the two unweighted Gaussians cross (used for reference / plotting)
% matches the PI's paper which does not apply mixture weights to the intersection

gaussPDF = @(x, mu, sig) (1 / (sig * sqrt(2*pi))) * exp(-(x - mu).^2 / (2 * sig^2));

intersectionEquation = @(x) gaussPDF(x, means(1), stdDevs(1)) - gaussPDF(x, means(2), stdDevs(2));
intersectionPoint = fzero(intersectionEquation, mean(means));

end


function [lowerBound, upperBound] = calculateMisclassificationBounds(means, stdDevs, thresh)
% Finds two asymmetric classification boundaries matching Miri et al. 2017 (Neuron) exactly.
% Equations from the paper (no mixture weights — assumes equal priors):

%   lowerBound b_L : mN(b_L) <= thresh
%     fraction of "narrow"-classified neurons that are actually wide
%     mN(b) = cdfW(b) / (cdfW(b) + cdfN(b))

%   upperBound b_U : mW(b_U) <= thresh
%     fraction of "wide"-classified neurons that are actually narrow
%     mW(b) = (1-cdfN(b)) / ((1-cdfN(b)) + (1-cdfW(b)))

% means(1), stdDevs(1) = narrow distribution (smaller mean)
% means(2), stdDevs(2) = wide distribution (larger mean)

muN = means(1); sigN = stdDevs(1);
muW = means(2); sigW = stdDevs(2);

cdfN = @(x) 0.5 * (1 + erf((x - muN) / (sigN * sqrt(2))));
cdfW = @(x) 0.5 * (1 + erf((x - muW) / (sigW * sqrt(2))));

% mN(b): fraction of narrow-classified that are actually wide
mN_eq = @(b) cdfW(b) / (cdfW(b) + cdfN(b)) - thresh;

% mW(b): fraction of wide-classified that are actually narrow
mW_eq = @(b) (1 - cdfN(b)) / ((1 - cdfN(b)) + (1 - cdfW(b))) - thresh;

x0 = mean(means);

bNarrowEq = fzero(mN_eq, x0);  % root of the narrow-side equation
bWideEq = fzero(mW_eq, x0);  % root of the wide-side equation

% fzero can converge to either side depending on equation shape, so don't
% assume which root is smaller — sort them explicitly to get the true
% lower (narrow cutoff) and upper (wide cutoff) bounds, guaranteeing a
% gap between them whenever one exists
lowerBound = min(bNarrowEq, bWideEq);
upperBound = max(bNarrowEq, bWideEq);

end
