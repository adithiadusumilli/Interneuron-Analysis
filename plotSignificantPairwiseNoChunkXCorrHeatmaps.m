function plotSignificantPairwiseNoChunkXCorrHeatmaps(combinedMatFile)
% heatmaps of significant pairwise NO-CHUNK peak correlations per session
% uses combined output:
%   all_peakCorrMat{sess}       (nInt x nPyr)
%   all_nullCorrMatShifts{sess} (nInt x nPyr x nShifts)

% significance: peakCorr > mean(null) + 3*std(null)
% nonsignificant entries set to NaN and rendered transparent

load(combinedMatFile, 'all_peakCorrMat','all_nullCorrMatShifts','baseDirs','all_nShifts');

numSessions = numel(all_peakCorrMat);

for sess = 1:numSessions
    peakCorrs = all_peakCorrMat{sess};
    nullXC    = all_nullCorrMatShifts{sess};

    if isempty(peakCorrs) || isempty(nullXC)
        fprintf('sess %d: empty data, skipping\n', sess);
        continue;
    end

    [nInt, nPyr, nS] = size(nullXC);
    if ~isequal(size(peakCorrs), [nInt nPyr])
        warning('sess %d: size mismatch (peakCorrs vs nullXC), skipping', sess);
        continue;
    end

    if exist('all_nShifts','var') && numel(all_nShifts) >= sess
        if all_nShifts(sess) ~= nS
            fprintf('sess %d: warning all_nShifts=%d but nullXC has %d\n', sess, all_nShifts(sess), nS);
        end
    end

    sigMat = nan(nInt, nPyr);
    sigPairs = 0;
    totalPairs = nInt * nPyr;

    for i = 1:nInt
        for j = 1:nPyr
            nullVals = squeeze(nullXC(i,j,:));
            if all(isnan(nullVals)), continue; end

            nullMean = mean(nullVals, 'omitnan');
            nullStd = std(nullVals,  'omitnan');

            a = peakCorrs(i,j);

            if ~isnan(a) && ~isnan(nullMean) && ~isnan(nullStd) && a > (nullMean + 3*nullStd)
                sigMat(i,j) = a;
                sigPairs = sigPairs + 1;
            end
        end
    end

    baseLabel = sprintf('sess %d', sess);
    if exist('baseDirs','var') && numel(baseDirs) >= sess
        baseLabel = sprintf('sess %d – %s', sess, baseDirs{sess});
    end

    figure('Name', sprintf('sess %d – NO-CHUNK significant peak corr', sess), 'Color', 'w');

    hImg = imagesc(sigMat);
    set(hImg, 'AlphaData', ~isnan(sigMat)); % hide nonsig pairs
    set(gca, 'Color', 'k'); % black background for NaNs

    axis xy;
    colormap(parula);
    c = colorbar;
    ylabel(c, 'peak correlation');

    xlabel('pyramidal neurons');
    ylabel('interneurons');

    title(sprintf('%s – NO-CHUNK significant pairs: %d / %d', baseLabel, sigPairs, totalPairs));
    set(gca, 'TickDir', 'out');
    box off;
end
end
