function plotSignificantPairwiseChunkedAllPlots_FDR_intPyr(fdrMatFile)
% makes heatmaps, histograms, and scatter plots for CHUNKED pairwise FDR results
% using only the interneuron x pyramidal block

% expects:
%   FDRresults.sessions{s}.peakCorrMat_all
%   FDRresults.sessions{s}.peakLagSecMat_all
%   FDRresults.sessions{s}.sigMaskFDR
%   FDRresults.sessions{s}.intIdx
%   FDRresults.sessions{s}.pyrIdx

S = load(fdrMatFile, 'FDRresults');
FDRresults = S.FDRresults;

numSessions = numel(FDRresults.sessions);

for sess = 1:numSessions
    D = FDRresults.sessions{sess};

    peakCorrs = D.peakCorrMat_all;
    peakLags = D.peakLagSecMat_all;
    sigMask = D.sigMaskFDR;

    if isempty(peakCorrs) || isempty(sigMask)
        fprintf('sess %d: empty data, skipping\n', sess);
        continue;
    end

    if isfield(D, 'intIdx') && ~isempty(D.intIdx)
        intIdx = D.intIdx(:);
    elseif isfield(D, 'neuronType') && ~isempty(D.neuronType)
        intIdx = find(D.neuronType == 1);
    else
        warning('sess %d: could not find intIdx, skipping', sess);
        continue;
    end

    if isfield(D, 'pyrIdx') && ~isempty(D.pyrIdx)
        pyrIdx = D.pyrIdx(:);
    elseif isfield(D, 'neuronType') && ~isempty(D.neuronType)
        pyrIdx = find(D.neuronType == 0);
    else
        warning('sess %d: could not find pyrIdx, skipping', sess);
        continue;
    end

    blockCorr = peakCorrs(intIdx, pyrIdx);
    blockLag = peakLags(intIdx, pyrIdx);
    blockSig = sigMask(intIdx, pyrIdx);

    sigCorrVec = blockCorr(blockSig);
    sigLagVec = blockLag(blockSig);

    allCorrVec = blockCorr(:);
    allLagVec = blockLag(:);

    sigCorrMat = blockCorr;
    sigCorrMat(~blockSig) = NaN;

    sigLagMat = blockLag;
    sigLagMat(~blockSig) = NaN;

    totalPairs = numel(allCorrVec);
    sigPairs = numel(sigCorrVec);

    animalID = D.animalID;
    if isempty(animalID)
        animalID = sprintf('sess %d', sess);
    end

    fprintf('\n=== %s ===\n', animalID);
    fprintf('int x pyr pairs: %d total | %d significant\n', totalPairs, sigPairs);

    % ---- lag bin edges ----
    if ~isempty(sigLagVec)
        minLagSec = floor(min(sigLagVec)*1000)/1000;
        maxLagSec = ceil(max(sigLagVec)*1000)/1000;
    elseif ~isempty(allLagVec)
        minLagSec = floor(min(allLagVec)*1000)/1000;
        maxLagSec = ceil(max(allLagVec)*1000)/1000;
    else
        minLagSec = -0.2;
        maxLagSec = 0.2;
    end
    lagBinEdgesSec = (minLagSec-0.001):0.001:(maxLagSec+0.001);

    % ==================== heatmap: peak correlation ====================
    figure('Name', sprintf('%s – CHUNKED int-pyr significant peak corr heatmap', animalID), 'Color', 'w');

    hImg = imagesc(sigCorrMat);
    set(hImg, 'AlphaData', ~isnan(sigCorrMat));
    set(gca, 'Color', 'k');

    axis xy;
    colormap(parula);
    c = colorbar;
    ylabel(c, 'peak correlation');

    xlabel('pyramidal neurons');
    ylabel('interneurons');
    title(sprintf('%s – CHUNKED int-pyr FDR-significant peak correlations: %d / %d', ...
        animalID, sigPairs, totalPairs));
    set(gca, 'TickDir', 'out');
    box off;

    % ==================== heatmap: peak lag ====================
    figure('Name', sprintf('%s – CHUNKED int-pyr significant peak lag heatmap', animalID), 'Color', 'w');

    hImg2 = imagesc(sigLagMat);
    set(hImg2, 'AlphaData', ~isnan(sigLagMat));
    set(gca, 'Color', 'k');

    axis xy;
    colormap(turbo);
    c2 = colorbar;
    ylabel(c2, 'peak lag (s)');

    xlabel('pyramidal neurons');
    ylabel('interneurons');
    title(sprintf('%s – CHUNKED int-pyr FDR-significant peak lags: %d / %d', ...
        animalID, sigPairs, totalPairs));
    set(gca, 'TickDir', 'out');
    box off;

    % ==================== histogram: significant peak lags ====================
    figure;
    histogram(sigLagVec, 'BinEdges', lagBinEdgesSec, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    xlabel('peak lag (s)');
    ylabel('count');
    title(sprintf('%s – CHUNKED int-pyr FDR-significant peak lags (n=%d / %d)', animalID, sigPairs, totalPairs));
    grid on;

    % ==================== histogram: significant peak correlations ====================
    figure;
    histogram(sigCorrVec, 50, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    xlabel('peak correlation');
    ylabel('count');
    title(sprintf('%s – CHUNKED int-pyr FDR-significant peak correlations (n=%d / %d)', animalID, sigPairs, totalPairs));
    grid on;

    % ==================== histogram: peak correlations > 0.2 ====================
    corrOverThresh = allCorrVec(allCorrVec > 0.2);
    figure;
    histogram(corrOverThresh, 50, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    xlabel('peak correlation (> 0.2)');
    ylabel('count');
    title(sprintf('%s – CHUNKED int-pyr peak correlations > 0.2 (n=%d / %d)', ...
        animalID, numel(corrOverThresh), totalPairs));
    grid on;

    % ==================== histogram: significant |lag| > 0.2 s ====================
    lagsOverThresh = sigLagVec(abs(sigLagVec) > 0.2);
    figure;
    histogram(lagsOverThresh, 'BinEdges', lagBinEdgesSec, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    xlabel('peak lag (s)');
    ylabel('count');
    title(sprintf('%s – CHUNKED int-pyr FDR-significant |peak lag| > 0.2 s (n=%d / %d sig)', ...
        animalID, numel(lagsOverThresh), numel(sigLagVec)));
    grid on;

    % ==================== scatterhist: significant pairs ====================
    if ~isempty(sigLagVec)
        figure;
        scatterhist(sigLagVec(:), sigCorrVec(:), 'Direction', 'out', 'Marker', '.');
        xlabel('peak lag (s)');
        ylabel('peak correlation');
        title(sprintf('%s – CHUNKED int-pyr FDR-significant pairs (n=%d / %d)', animalID, sigPairs, totalPairs));
    end

    % ==================== corr >= 0.2 analyses ====================
    corrMask = (allCorrVec >= 0.2) & ~isnan(allCorrVec) & ~isnan(allLagVec);
    corrOverThresh = allCorrVec(corrMask);
    lagsForCorrOverThresh = allLagVec(corrMask);

    figure;
    histogram(corrOverThresh, 50, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    xlabel('peak correlation (>= 0.2)');
    ylabel('count');
    title(sprintf('%s – CHUNKED int-pyr peak correlations >= 0.2 (n=%d / %d)', ...
        animalID, numel(corrOverThresh), totalPairs));
    grid on;

    if ~isempty(lagsForCorrOverThresh)
        minLagSec2 = floor(min(lagsForCorrOverThresh)*1000)/1000;
        maxLagSec2 = ceil(max(lagsForCorrOverThresh)*1000)/1000;
        lagBinEdgesSec2 = (minLagSec2-0.001):0.001:(maxLagSec2+0.001);

        figure;
        histogram(lagsForCorrOverThresh, 'BinEdges', lagBinEdgesSec2, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
        xlabel('peak lag (s) for pairs with corr >= 0.2');
        ylabel('count');
        title(sprintf('%s – CHUNKED int-pyr peak lags where peak corr >= 0.2 (n=%d)', ...
            animalID, numel(lagsForCorrOverThresh)));
        grid on;

        figure;
        scatter(lagsForCorrOverThresh(:), corrOverThresh(:), '.');
        xlabel('peak lag (s)');
        ylabel('peak correlation');
        title(sprintf('%s – CHUNKED int-pyr lag vs corr for corr >= 0.2 (n=%d)', animalID, numel(corrOverThresh)));
        grid on;
    end
end
end
