function plotSignificantPairwiseNoChunkXCorrsPerAnimal(combinedMatFile)
% plot significant pairwise NO-CHUNK xcorrs per animal (from quest runs) using combined output:
%   all_peakLagMat{sess}      (nInt x nPyr)   real peak lag (sec)
%   all_peakCorrMat{sess}     (nInt x nPyr)   real peak corr
%   all_nullCorrShifts{sess}  (nInt x nPyr x nShifts)   null corr (0-lag) per shift

% significance: peakCorr > mean(null) + 3*std(null)

% also plots:
%   - peak corr > 0.2 (no significance)
%   - significant |lag| > 0.2s
%   - scatterhist of significant pairs

load(combinedMatFile, 'all_peakLagMat','all_peakCorrMat','all_nullCorrShifts','baseDirs','nShifts');

numSessions = numel(all_peakCorrMat);

for sess = 1:numSessions
    peakCorrs = all_peakCorrMat{sess};
    peakLags = all_peakLagMat{sess};
    nullXC = all_nullCorrShifts{sess}; % nInt x nPyr x nShifts

    if isempty(peakCorrs) || isempty(nullXC)
        fprintf('sess %d: empty data, skipping\n', sess);
        continue;
    end

    if ndims(nullXC) ~= 3
        warning('sess %d: nullXC is not 3D, skipping', sess);
        continue;
    end

    [nInt, nPyr, nS] = size(nullXC);
    if exist('nShifts','var') && ~isempty(nShifts) && nS ~= nShifts
        fprintf('sess %d: warning nShifts in file=%d but nullXC has %d\n', sess, nShifts, nS);
    end

    if ~isequal(size(peakCorrs), [nInt nPyr]) || ~isequal(size(peakLags), [nInt nPyr])
        warning('sess %d: size mismatch between peak mats and nullXC, skipping', sess);
        continue;
    end

    allCorrVec = peakCorrs(:);
    allLagVec  = peakLags(:);

    sigMask = false(nInt, nPyr);

    % compute significance per pair
    for i = 1:nInt
        for j = 1:nPyr
            nullVals = squeeze(nullXC(i,j,:));
            nullMean = mean(nullVals, 'omitnan');
            nullStd = std(nullVals, 'omitnan');

            a = peakCorrs(i,j);

            if ~isnan(a) && ~isnan(nullMean) && ~isnan(nullStd) && a > (nullMean + 3*nullStd)
                sigMask(i,j) = true;
            end
        end
    end

    sigCorrVec = peakCorrs(sigMask);
    sigLagVec = peakLags(sigMask);

    totalPairs = numel(allCorrVec);
    sigPairs = numel(sigCorrVec);

    fprintf('\n=== sess %d ===\n', sess);
    if exist('baseDirs','var') && numel(baseDirs) >= sess
        fprintf('baseDir: %s\n', baseDirs{sess});
    end
    fprintf('pairs: %d total | %d significant\n', totalPairs, sigPairs);

    % ---- 1 ms bin edges for lag histogram (sec) ----
    if ~isempty(sigLagVec)
        minLagSec = floor(min(sigLagVec)*1000)/1000;
        maxLagSec = ceil(max(sigLagVec)*1000)/1000;
    elseif ~isempty(allLagVec)
        minLagSec = floor(min(allLagVec)*1000)/1000;
        maxLagSec = ceil(max(allLagVec)*1000)/1000;
    else
        minLagSec = -0.2; maxLagSec = 0.2;
    end
    lagBinEdgesSec = (minLagSec-0.001):0.001:(maxLagSec+0.001);

    % ---- plot: significant peak lags ----
    figure;
    histogram(sigLagVec, 'BinEdges', lagBinEdgesSec, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    xlabel('peak lag (s)');
    ylabel('count');
    title(sprintf('sess %d – NO-CHUNK significant peak lags (n=%d / %d)', sess, sigPairs, totalPairs));
    grid on;

    % ---- plot: significant peak correlations ----
    figure;
    histogram(sigCorrVec, 50, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    xlabel('peak correlation');
    ylabel('count');
    title(sprintf('sess %d – NO-CHUNK significant peak correlations (n=%d / %d)', sess, sigPairs, totalPairs));
    grid on;

    % ---- plot: peak correlations > 0.2 (no significance) ----
    corrOverThresh = allCorrVec(allCorrVec > 0.2);
    figure;
    histogram(corrOverThresh, 50, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    xlabel('peak correlation (> 0.2)');
    ylabel('count');
    title(sprintf('sess %d – NO-CHUNK peak correlations > 0.2 (n=%d / %d)', sess, numel(corrOverThresh), totalPairs));
    grid on;

    % ---- plot: significant |lag| > 0.2 s ----
    lagsOverThresh = sigLagVec(abs(sigLagVec) > 0.2);
    figure;
    histogram(lagsOverThresh, 'BinEdges', lagBinEdgesSec, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    xlabel('peak lag (s)');
    ylabel('count');
    title(sprintf('sess %d – NO-CHUNK significant |peak lag| > 0.2 s (n=%d / %d sig)', ...
        sess, numel(lagsOverThresh), numel(sigLagVec)));
    grid on;

    % ---- scatterhist: significant pairs ----
    if ~isempty(sigLagVec)
        figure;
        scatterhist(sigLagVec(:), sigCorrVec(:), 'Direction', 'out', 'Marker', '.');
        xlabel('peak lag (s)');
        ylabel('peak correlation');
        title(sprintf('sess %d – NO-CHUNK significant pairs (n=%d / %d)', sess, sigPairs, totalPairs));
    end
end
end
