function classifyCortexAndStriatum(neuronDataStructFiles, firingRatesFiles, outFile)
% performs gmm-based spike-width classification for BOTH cortex and striatum and saves combined results in one .mat file

% inputs:
%   neuronDataStructFiles : cell array of full paths to neuronDataStruct.mat
%   firingRatesFiles      : cell array of full paths to NeuralFiringRates*.mat
%   outFile               : (optional) name for output .mat  [default = 'AA_classifications.mat']

% output (saved):
%   classifications{fileIdx,1}   → cortex labels  (0=pyramidal, 1=interneuron, NaN=other)
%   classifications{fileIdx,2}   → striatum labels (same coding)

if nargin < 3
    outFile = 'AA_classifications.mat';
end

conslidatedDataFoler = 'X:\David\AnalysesData';

numFiles        = numel(neuronDataStructFiles);
classifications = cell(numFiles,2);      % col-1 cortex, col-2 striatum
fileWidths      = cell(numFiles,2);      % store per-file widths for later

% PART 1: determine spike-widths for cortex and striatum 

allWidths = {[],[]};   % {1}=cortex, {2}=striatum

for fileIndex = 1:numFiles
    load(neuronDataStructFiles{fileIndex},'neuronDataStruct');
    load(firingRatesFiles{fileIndex},'cortexInds','striatumInds');

    % -------- cortex
    cortexWidths = [];
    for i = 1:numel(cortexInds)
        waveform    = neuronDataStruct(cortexInds(i)).waveforms;
        biggestChan = neuronDataStruct(cortexInds(i)).biggestChan;
        ap          = waveform(:,biggestChan);
        [~,mx] = max(ap); [~,mn] = min(ap);
        cortexWidths = [cortexWidths , abs(mx-mn)/30000];  % sec
    end
    fileWidths{fileIndex,1} = cortexWidths;
    allWidths{1}            = [allWidths{1} , cortexWidths];

    % -------- striatum
    striatWidths = [];
    for i = 1:numel(striatumInds)
        waveform    = neuronDataStruct(striatumInds(i)).waveforms;
        biggestChan = neuronDataStruct(striatumInds(i)).biggestChan;
        ap          = waveform(:,biggestChan);
        [~,mx] = max(ap); [~,mn] = min(ap);
        striatWidths = [striatWidths , abs(mx-mn)/30000];  % sec
    end
    fileWidths{fileIndex,2} = striatWidths;
    allWidths{2}            = [allWidths{2} , striatWidths];
end

% PARTS 2 & 3: fit GMMs & plot spike widths for each region per animal
numComponents = 2;
intersectionPointCortex = nan(1, numFiles);
intersectionPointStriat = nan(1, numFiles);
for fileIndex = 1:numFiles
    
    % -----Cortex-----
    cortexWidths = fileWidths{fileIndex, 1};
    gm = fitgmdist(cortexWidths', numComponents);  % Use per-animal widths
    means = gm.mu;
    stdDevs = sqrt(gm.Sigma);
    weights = gm.ComponentProportion;
    % Sort components: narrow (INT) first, wide (PYR) second
    [means, order] = sort(means, 'ascend');
    stdDevs = stdDevs(order);
    weights = weights(order);
    % Calculate per-animal intersection BEFORE plotting
    intersectionPointCortex(fileIndex) = calculateIntersectionPoint(means, stdDevs);
    % Plot
    figure;
    h = histogram(cortexWidths, 'BinWidth', 1/20000, 'EdgeColor', 'black', 'FaceColor', 'blue'); hold on;
    x = linspace(min(cortexWidths), max(cortexWidths), 1000);
    yInt = pdf('Normal', x, means(1), stdDevs(1)) * numel(cortexWidths) * h.BinWidth * weights(1);
    yPyr = pdf('Normal', x, means(2), stdDevs(2)) * numel(cortexWidths) * h.BinWidth * weights(2);
    plot(x, yPyr, 'r', 'LineWidth', 2);   % Pyramidal = wide
    plot(x, yInt, 'g', 'LineWidth', 2);   % Interneuron = narrow
    xline(intersectionPointCortex(fileIndex), 'k--');  % Use correct line
    title(sprintf('Animal %d: Cortex Spike Widths with Fitted GMMs', fileIndex));
    legend({'Widths','Pyramidal','Interneuron'});
    box off;
    % -----Striatum-----
    striatWidths = fileWidths{fileIndex, 2};
    gm = fitgmdist(striatWidths', numComponents);  % Use per-animal widths
    means = gm.mu;
    stdDevs = sqrt(gm.Sigma);
    weights = gm.ComponentProportion;
    % Sort components
    [means, order] = sort(means, 'ascend');
    stdDevs = stdDevs(order);
    weights = weights(order);
    % Calculate intersection point
    intersectionPointStriat(fileIndex) = calculateIntersectionPoint(means, stdDevs);
    % Plot
    figure;
    h = histogram(striatWidths, 'BinWidth', 1/20000, 'EdgeColor', 'black', 'FaceColor', 'magenta'); hold on;
    x = linspace(min(striatWidths), max(striatWidths), 1000);
    yInt = pdf('Normal', x, means(1), stdDevs(1)) * numel(striatWidths) * h.BinWidth * weights(1);
    yPyr = pdf('Normal', x, means(2), stdDevs(2)) * numel(striatWidths) * h.BinWidth * weights(2);
    plot(x, yPyr, 'r', 'LineWidth', 2);   % Pyramidal = wide
    plot(x, yInt, 'g', 'LineWidth', 2);   % Interneuron = narrow
    xline(intersectionPointStriat(fileIndex), 'k--');  % Correct line
    title(sprintf('Animal %d: Striatum Spike Widths with Fitted GMMs', fileIndex));
    legend({'Widths','Pyramidal','Interneuron'});
    box off;
end

% PART 4. assign labels per file

for fileIndex = 1:numFiles
    load(neuronDataStructFiles{fileIndex},'neuronDataStruct');
    load(firingRatesFiles{fileIndex},'cortexInds','striatumInds');

    % initialize arrays
    cortexLabels   = nan(1,numel(neuronDataStruct));
    striatLabels   = nan(1,numel(neuronDataStruct));

    % cortex
    for i = 1:numel(cortexInds)
        if fileWidths{fileIndex,1}(i) < intersectionPointCortex
            cortexLabels(cortexInds(i)) = 1;  % interneuron
        else
            cortexLabels(cortexInds(i)) = 0;  % pyramidal
        end
    end

    % striatum
    for i = 1:numel(striatumInds)
        if fileWidths{fileIndex,2}(i) < intersectionPointStriat
            striatLabels(striatumInds(i)) = 1; % interneuron-like
        else
            striatLabels(striatumInds(i)) = 0; % pyramidal-like (MSN)
        end
    end

    classifications{fileIndex,1} = cortexLabels;
    classifications{fileIndex,2} = striatLabels;
end

% PART 5. save & summary

save(fullfile(conslidatedDataFoler,outFile),'classifications');
fprintf('saved combined cortex & striatum classifications to %s\n',outFile);

for f = 1:numFiles
    fprintf('file %d: cortex  pyr=%d  int=%d | striatum  pyr=%d  int=%d\n',...
        f, sum(classifications{f,1}==0,'omitnan'), sum(classifications{f,1}==1,'omitnan'), ...
           sum(classifications{f,2}==0,'omitnan'), sum(classifications{f,2}==1,'omitnan'));
end
end

% EXTRA: helper function for calculating intersection point

function intersectionPoint = calculateIntersectionPoint(means,stdDevs)
gaussPDF = @(x,mu,sig) (1/(sig*sqrt(2*pi))) * exp(-(x-mu).^2/(2*sig^2));
intersectionEquation = @(x) gaussPDF(x,means(1),stdDevs(1)) - gaussPDF(x,means(2),stdDevs(2));
intersectionPoint = fzero(intersectionEquation,mean(means));
end
