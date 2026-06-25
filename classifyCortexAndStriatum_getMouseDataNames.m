function classifyCortexAndStriatum_getMouseDataNames(mouseIDs, baseSessionNames, probeRegions, outFile)
% classifies cortex and striatum neurons using waveform width + GMM
% this version uses David's getMouseDataNames cuz new animals have some diff file naming

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

%% make sure MATLAB can find getMouseDataNames.m
addpath('C:\Users\mirilab\Documents\GlobusTransfer');

consolidatedDataFolder = 'X:\David\AnalysesData';

numFiles = numel(mouseIDs);

if numel(baseSessionNames) ~= numFiles || numel(probeRegions) ~= numFiles
    error('mouseIDs, baseSessionNames, and probeRegions must all have the same length');
end

classifications = cell(numFiles,2); % col 1 = cortex, col 2 = striatum
fileWidths = cell(numFiles,2); % col 1 = cortex widths, col 2 = striatum widths

neuronDataStructFiles = cell(numFiles,1);
firingRatesFiles = cell(numFiles,1);

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
    fprintf('firing rates:      %s\n', firingRatesFiles{fileIndex});
end

%% ---------------- PART 1: calculate spike widths ----------------
% width = time between waveform peak and trough on the biggest channel

for fileIndex = 1:numFiles

    load(neuronDataStructFiles{fileIndex}, 'neuronDataStruct');
    load(firingRatesFiles{fileIndex}, 'cortexInds', 'striatumInds');

    %% cortex widths
    cortexWidths = nan(1,numel(cortexInds));

    for i = 1:numel(cortexInds)
        waveform = neuronDataStruct(cortexInds(i)).waveforms;
        biggestChan = neuronDataStruct(cortexInds(i)).biggestChan;

        ap = waveform(:,biggestChan);

        [~, mx] = max(ap);
        [~, mn] = min(ap);

        cortexWidths(i) = abs(mx - mn) / 30000;  % seconds
    end

    fileWidths{fileIndex,1} = cortexWidths;

    %% striatum widths
    striatWidths = nan(1,numel(striatumInds));

    for i = 1:numel(striatumInds)
        waveform = neuronDataStruct(striatumInds(i)).waveforms;
        biggestChan = neuronDataStruct(striatumInds(i)).biggestChan;

        ap = waveform(:,biggestChan);

        [~, mx] = max(ap);
        [~, mn] = min(ap);

        striatWidths(i) = abs(mx - mn) / 30000;  % seconds
    end

    fileWidths{fileIndex,2} = striatWidths;
end

%% ---------------- PART 2: fit GMMs and find per-animal cutoffs ----------------
numComponents = 2;

intersectionPointCortex = nan(1,numFiles);
intersectionPointStriat = nan(1,numFiles);

for fileIndex = 1:numFiles

    animalID = mouseIDs{fileIndex};

    %% cortex GMM
    cortexWidths = fileWidths{fileIndex,1};
    cortexWidths = cortexWidths(~isnan(cortexWidths) & isfinite(cortexWidths));

    gm = fitgmdist(cortexWidths', numComponents);

    means = gm.mu;
    stdDevs = squeeze(sqrt(gm.Sigma));
    weights = gm.ComponentProportion;

    [means, order] = sort(means, 'ascend');  % narrow first, wide second
    stdDevs = stdDevs(order);
    weights = weights(order);

    intersectionPointCortex(fileIndex) = calculateIntersectionPoint(means, stdDevs);

    figure;
    h = histogram(cortexWidths, 'BinWidth', 1/20000, ...
        'EdgeColor', 'black', 'FaceColor', 'blue');
    hold on;

    x = linspace(min(cortexWidths), max(cortexWidths), 1000);

    yInt = pdf('Normal', x, means(1), stdDevs(1)) * numel(cortexWidths) * h.BinWidth * weights(1);
    yPyr = pdf('Normal', x, means(2), stdDevs(2)) * numel(cortexWidths) * h.BinWidth * weights(2);

    plot(x, yPyr, 'r', 'LineWidth', 2);
    plot(x, yInt, 'g', 'LineWidth', 2);
    xline(intersectionPointCortex(fileIndex), 'k--', 'LineWidth', 1.5);

    title(sprintf('%s Cortex Spike Widths with GMM Cutoff', animalID));
    legend({'Widths','Pyramidal','Interneuron','Cutoff'});
    box off;

    %% striatum GMM
    striatWidths = fileWidths{fileIndex,2};
    striatWidths = striatWidths(~isnan(striatWidths) & isfinite(striatWidths));

    gm = fitgmdist(striatWidths', numComponents);

    means = gm.mu;
    stdDevs = squeeze(sqrt(gm.Sigma));
    weights = gm.ComponentProportion;

    [means, order] = sort(means, 'ascend');  % narrow first, wide second
    stdDevs = stdDevs(order);
    weights = weights(order);

    intersectionPointStriat(fileIndex) = calculateIntersectionPoint(means, stdDevs);

    figure;
    h = histogram(striatWidths, 'BinWidth', 1/20000, ...
        'EdgeColor', 'black', 'FaceColor', 'magenta');
    hold on;

    x = linspace(min(striatWidths), max(striatWidths), 1000);

    yInt = pdf('Normal', x, means(1), stdDevs(1)) * numel(striatWidths) * h.BinWidth * weights(1);
    yPyr = pdf('Normal', x, means(2), stdDevs(2)) * numel(striatWidths) * h.BinWidth * weights(2);

    plot(x, yPyr, 'r', 'LineWidth', 2);
    plot(x, yInt, 'g', 'LineWidth', 2);
    xline(intersectionPointStriat(fileIndex), 'k--', 'LineWidth', 1.5);

    title(sprintf('%s Striatum Spike Widths with GMM Cutoff', animalID));
    legend({'Widths','Wide','Narrow','Cutoff'});
    box off;
end

%% ---------------- PART 3: assign labels using each animal's own cutoff ----------------
% label convention:
%   0 = pyramidal / wide waveform
%   1 = interneuron / narrow waveform
%   NaN = not in that region

for fileIndex = 1:numFiles

    load(neuronDataStructFiles{fileIndex}, 'neuronDataStruct');
    load(firingRatesFiles{fileIndex}, 'cortexInds', 'striatumInds');

    cortexLabels = nan(1,numel(neuronDataStruct));
    striatLabels = nan(1,numel(neuronDataStruct));

    %% cortex labels
    for i = 1:numel(cortexInds)

        if fileWidths{fileIndex,1}(i) < intersectionPointCortex(fileIndex)
            cortexLabels(cortexInds(i)) = 1;  % narrow = interneuron
        else
            cortexLabels(cortexInds(i)) = 0;  % wide = pyramidal
        end
    end

    %% striatum labels
    for i = 1:numel(striatumInds)

        if fileWidths{fileIndex,2}(i) < intersectionPointStriat(fileIndex)
            striatLabels(striatumInds(i)) = 1;  % narrow
        else
            striatLabels(striatumInds(i)) = 0;  % wide
        end
    end

    classifications{fileIndex,1} = cortexLabels;
    classifications{fileIndex,2} = striatLabels;
end

%% ---------------- PART 4: save everything ----------------
outPath = fullfile(consolidatedDataFolder, outFile);

save(outPath, ...
    'classifications', ...
    'mouseIDs', ...
    'baseSessionNames', ...
    'probeRegions', ...
    'neuronDataStructFiles', ...
    'firingRatesFiles', ...
    'fileWidths', ...
    'intersectionPointCortex', ...
    'intersectionPointStriat', ...
    '-v7.3');

fprintf('\nsaved combined cortex & striatum classifications to:\n%s\n', outPath);

for f = 1:numFiles
    fprintf('%s: cortex pyr=%d int=%d | striatum wide=%d narrow=%d\n', ...
        mouseIDs{f}, ...
        sum(classifications{f,1}==0,'omitnan'), ...
        sum(classifications{f,1}==1,'omitnan'), ...
        sum(classifications{f,2}==0,'omitnan'), ...
        sum(classifications{f,2}==1,'omitnan'));
end

end

%% ========================================================================
%% helper
%% ========================================================================
function intersectionPoint = calculateIntersectionPoint(means,stdDevs)

gaussPDF = @(x,mu,sig) ...
    (1/(sig*sqrt(2*pi))) * exp(-(x-mu).^2/(2*sig^2));

intersectionEquation = @(x) ...
    gaussPDF(x,means(1),stdDevs(1)) - gaussPDF(x,means(2),stdDevs(2));

intersectionPoint = fzero(intersectionEquation, mean(means));

end
