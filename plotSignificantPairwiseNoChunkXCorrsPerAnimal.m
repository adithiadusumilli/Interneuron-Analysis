function plotSignificantPairwiseNoChunkXCorrsPerAnimal(combinedMatFile)
% plot significant pairwise NO-CHUNK xcorrs per animal using combined output:
%   all_peakLagMat{sess}        (nInt x nPyr)  in seconds
%   all_peakCorrMat{sess}       (nInt x nPyr)
%   all_nullCorrMatShifts{sess} (nInt x nPyr x nShifts)

% significance: peakCorr > mean(null) + 3*std(null)

% also plots:
%   - peak corr > 0.2 (no significance)
%   - significant |lag| > 0.2s
%   - scatterhist of significant pairs

load(combinedMatFile, ...
    'all_peakLagMat','all_peakCorrMat','all_nullCorrMatShifts', ...
    'baseDirs','all_nShifts');

numSessions = numel(all_peakCorrMat);

for sess = 1:numSessions
    peakCorrs = all_peakCorrMat{sess};
    peakLags = all_peakLagMat{sess};
    nullXC = all_nullCorrMatShifts{sess}; % nInt x nPyr x nShifts

    if isempty(peakCorrs) || isempty(nullXC)
        fprintf('sess %d: empty data, skipping\n', sess);
        continue;
    end

    if ndims(nullXC) ~= 3
        warning('sess %d: nullXC is not 3D, skipping', sess);
        continue;
    end

    [nInt, nPyr, nS] = size(nullXC);

    % nShifts for this session (if present)
    if exist('all_nShifts','var') && numel(all_nShifts) >= sess && ~isempty(all_nShifts(sess))
        nShifts_expected = all_nShifts(sess);
        if nShifts_expected ~= nS
            fprintf('sess %d: warning all_nShifts=%d but nullXC has %d\n', sess, nShifts_expected, nS);
        end
    end

    if ~isequal(size(peakCorrs), [nInt nPyr]) || ~isequal(size(peakLags), [nInt nPyr])
        warning('sess %d: size mismatch between peak mats and nullXC, skipping', sess);
        continue;
    end

    % compute significance per pair
    sigMask = false(nInt, nPyr);

    for i = 1:nInt
        for j = 1:nPyr
            nullVals = squeeze(nullXC(i,j,:));
            nullMean = mean(nullVals, 'omitnan');
            nullStd = std(nullVals,  'omitnan');

            a = peakCorrs(i,j);

            if ~isnan(a) && ~isnan(nullMean) && ~isnan(nullStd) && a > (nullMean + 3*nullStd)
                sigMask(i,j) = true;
            end
        end
    end

    sigCorrVec = peakCorrs(sigMask);
    sigLagVec = peakLags(sigMask);

    allCorrVec = peakCorrs(:);
    allLagVec = peakLags(:);

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
    title(sprintf('sess %d – NO-CHUNK peak correlations > 0.2 (n=%d / %d)', ...
        sess, numel(corrOverThresh), totalPairs));
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

    % ---- plot: peak correlations >= 0.2 (no significance) ----
corrMask = (allCorrVec >= 0.2) & ~isnan(allCorrVec) & ~isnan(allLagVec);
corrOverThresh = allCorrVec(corrMask);
lagsForCorrOverThresh = allLagVec(corrMask);

figure;
histogram(corrOverThresh, 50, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
xlabel('peak correlation (>= 0.2)');
ylabel('count');
title(sprintf('sess %d – NO-CHUNK peak correlations >= 0.2 (n=%d / %d)', ...
    sess, numel(corrOverThresh), totalPairs));
grid on;

% new plot the peak LAGS corresponding to corr >= 0.2
if ~isempty(lagsForCorrOverThresh)
    % 1 ms bin edges (sec) for these lags
    minLagSec2 = floor(min(lagsForCorrOverThresh)*1000)/1000;
    maxLagSec2 = ceil(max(lagsForCorrOverThresh)*1000)/1000;
    lagBinEdgesSec2 = (minLagSec2-0.001):0.001:(maxLagSec2+0.001);

    figure;
    histogram(lagsForCorrOverThresh, 'BinEdges', lagBinEdgesSec2, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    xlabel('peak lag (s) for pairs with corr >= 0.2');
    ylabel('count');
    title(sprintf('sess %d – NO-CHUNK peak lags where peak corr >= 0.2 (n=%d)', ...
        sess, numel(lagsForCorrOverThresh)));
    grid on;

    % scatter lag vs corr for corr>=0.2
    figure;
    scatter(lagsForCorrOverThresh(:), corrOverThresh(:), '.');
    xlabel('peak lag (s)');
    ylabel('peak correlation');
    title(sprintf('sess %d – NO-CHUNK lag vs corr for corr >= 0.2 (n=%d)', sess, numel(corrOverThresh)));
    grid on;
end

end
end
