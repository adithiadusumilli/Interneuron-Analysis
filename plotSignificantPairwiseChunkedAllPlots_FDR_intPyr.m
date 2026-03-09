function plotSignificantPairwiseChunkedAllPlots_FDR_intPyr(fdrMatFile)
% makes heatmaps, histograms, and scatter plots for CHUNKED pairwise FDR results
% using only the interneuron x pyramidal block

% expects:
%   FDRresults.sessions{s}.peakCorrMat_all
%   FDRresults.sessions{s}.peakLagSecMat_all
%   FDRresults.sessions{s}.sigMaskFDR
%   FDRresults.sessions{s}.nInt_ref
%   FDRresults.sessions{s}.nPyr_ref

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

    if ~isfield(D, 'nInt_ref') || ~isfield(D, 'nPyr_ref') || isempty(D.nInt_ref) || isempty(D.nPyr_ref)
        warning('sess %d: missing nInt_ref / nPyr_ref, skipping', sess);
        continue;
    end

    % use matrix ordering directly: interneurons first, pyramidal second
    intRows = 1:D.nInt_ref;
    pyrCols = D.nInt_ref + (1:D.nPyr_ref);

    if max(pyrCols) > size(peakCorrs,2)
        warning('sess %d: int/pyr block exceeds matrix dimensions, skipping', sess);
        continue;
    end

    blockCorr = peakCorrs(intRows, pyrCols);
    blockLag = peakLags(intRows, pyrCols);
    blockSig = sigMask(intRows, pyrCols);

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
    fprintf('%s: heatmap non-nan entries = %d\n', animalID, nnz(~isnan(sigCorrMat)));

    % ---- robust lag bin edges ----
    lagSource = sigLagVec(~isnan(sigLagVec) & isfinite(sigLagVec));
    if isempty(lagSource)
        lagSource = allLagVec(~isnan(allLagVec) & isfinite(allLagVec));
    end

    if isempty(lagSource)
        lagBinEdgesSec = -0.201:0.001:0.201;
    else
        minLagSec = floor(min(lagSource)*1000)/1000;
        maxLagSec = ceil(max(lagSource)*1000)/1000;

        if minLagSec == maxLagSec
            minLagSec = minLagSec - 0.001;
            maxLagSec = maxLagSec + 0.001;
        end

        lagBinEdgesSec = minLagSec:0.001:maxLagSec;

        if numel(lagBinEdgesSec) < 2
            lagBinEdgesSec = [minLagSec-0.001, maxLagSec+0.001];
        end
    end

    % ==================== reference heatmap: all int x pyr peak correlation ====================
    %figure('Name', sprintf('%s – CHUNKED int-pyr all peak corr heatmap', animalID), 'Color', 'w');
    %imagesc(blockCorr);
    %axis xy;
    %colormap(parula);
    %c = colorbar;
    %ylabel(c, 'peak correlation');
    %xlabel('pyramidal neurons');
    %ylabel('interneurons');
    %title(sprintf('%s – CHUNKED int-pyr all peak correlations', animalID));
    %set(gca, 'TickDir', 'out');
    5box off;

    % ==================== significant heatmap: peak correlation ====================
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

    % ==================== reference heatmap: all int x pyr peak lag ====================
    %figure('Name', sprintf('%s – CHUNKED int-pyr all peak lag heatmap', animalID), 'Color', 'w');
    %imagesc(blockLag);
    %axis xy;
    %colormap(turbo);
    %c2 = colorbar;
    %ylabel(c2, 'peak lag (s)');
    %xlabel('pyramidal neurons');
    %ylabel('interneurons');
    %title(sprintf('%s – CHUNKED int-pyr all peak lags', animalID));
    %set(gca, 'TickDir', 'out');
    %box off;

    % ==================== significant heatmap: peak lag ====================
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

    if sigPairs == 0
        fprintf('%s: no FDR-significant int x pyr pairs, so significant histograms/scatters were skipped.\n', animalID);
    else
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

        % ==================== scatterhist: significant pairs ====================
        figure;
        scatterhist(sigLagVec(:), sigCorrVec(:), 'Direction', 'out', 'Marker', '.');
        xlabel('peak lag (s)');
        ylabel('peak correlation');
        title(sprintf('%s – CHUNKED int-pyr FDR-significant pairs (n=%d / %d)', animalID, sigPairs, totalPairs));
    end

    % ==================== histogram: peak correlations >= 0.2 ====================
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
        lagSource2 = lagsForCorrOverThresh(~isnan(lagsForCorrOverThresh) & isfinite(lagsForCorrOverThresh));

        if isempty(lagSource2)
            lagBinEdgesSec2 = -0.201:0.001:0.201;
        else
            minLagSec2 = floor(min(lagSource2)*1000)/1000;
            maxLagSec2 = ceil(max(lagSource2)*1000)/1000;

            if minLagSec2 == maxLagSec2
                minLagSec2 = minLagSec2 - 0.001;
                maxLagSec2 = maxLagSec2 + 0.001;
            end

            lagBinEdgesSec2 = minLagSec2:0.001:maxLagSec2;

            if numel(lagBinEdgesSec2) < 2
                lagBinEdgesSec2 = [minLagSec2-0.001, maxLagSec2+0.001];
            end
        end

        % ==================== histogram: peak lags where corr >= 0.2 ====================
        figure;
        histogram(lagsForCorrOverThresh, 'BinEdges', lagBinEdgesSec2, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
        xlabel('peak lag (s) for pairs with corr >= 0.2');
        ylabel('count');
        title(sprintf('%s – CHUNKED int-pyr peak lags where peak corr >= 0.2 (n=%d)', ...
            animalID, numel(lagsForCorrOverThresh)));
        grid on;

        % ==================== scatter: lag vs corr where corr >= 0.2 ====================
        figure;
        scatter(lagsForCorrOverThresh(:), corrOverThresh(:), '.');
        xlabel('peak lag (s)');
        ylabel('peak correlation');
        title(sprintf('%s – CHUNKED int-pyr lag vs corr for corr >= 0.2 (n=%d)', animalID, numel(corrOverThresh)));
        grid on;
    end
end
end
