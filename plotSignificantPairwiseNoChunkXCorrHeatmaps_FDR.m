function plotSignificantPairwiseNoChunkXCorrHeatmaps_FDR(fdrMatFile)
% heatmaps of significant pairwise NO-CHUNK peak correlations and peak lags per session
% using FDR results

% expects:
%   FDRresults.sessions{s}.peakCorrMatAll
%   FDRresults.sessions{s}.peakLagMatAll
%   FDRresults.sessions{s}.sigMaskFDR

% nonsignificant entries are set to NaN and rendered transparent

S = load(fdrMatFile, 'FDRresults');
FDRresults = S.FDRresults;

numSessions = numel(FDRresults.sessions);

for sess = 1:numSessions
    D = FDRresults.sessions{sess};

    peakCorrs = D.peakCorrMatAll;
    peakLags = D.peakLagMatAll;
    sigMask = D.sigMaskFDR;

    if isempty(peakCorrs) || isempty(sigMask)
        fprintf('sess %d: empty data, skipping\n', sess);
        continue;
    end

    if ~isequal(size(peakCorrs), size(sigMask)) || ~isequal(size(peakLags), size(sigMask))
        warning('sess %d: size mismatch, skipping', sess);
        continue;
    end

    sigCorrMat = peakCorrs;
    sigCorrMat(~sigMask) = NaN;

    sigLagMat = peakLags;
    sigLagMat(~sigMask) = NaN;

    upperPairs = nnz(triu(true(size(sigMask)),1));
    sigPairs = nnz(triu(sigMask,1));

    animalID = D.animalID;
    if isempty(animalID)
        animalID = sprintf('sess %d', sess);
    end

    % ---- heatmap: peak correlation ----
    figure('Name', sprintf('%s – NO-CHUNK significant peak corr heatmap', animalID), 'Color', 'w');

    hImg = imagesc(sigCorrMat);
    set(hImg, 'AlphaData', ~isnan(sigCorrMat));
    set(gca, 'Color', 'k');

    axis xy;
    colormap(parula);
    c = colorbar;
    ylabel(c, 'peak correlation');

    xlabel('neuron j');
    ylabel('neuron i');
    title(sprintf('%s – NO-CHUNK FDR-significant peak correlations: %d / %d', ...
        animalID, sigPairs, upperPairs));
    set(gca, 'TickDir', 'out');
    box off;

    % ---- heatmap: peak lag ----
    figure('Name', sprintf('%s – NO-CHUNK significant peak lag heatmap', animalID), 'Color', 'w');

    hImg2 = imagesc(sigLagMat);
    set(hImg2, 'AlphaData', ~isnan(sigLagMat));
    set(gca, 'Color', 'k');

    axis xy;
    colormap(turbo);
    c2 = colorbar;
    ylabel(c2, 'peak lag (s)');

    xlabel('neuron j');
    ylabel('neuron i');
    title(sprintf('%s – NO-CHUNK FDR-significant peak lags: %d / %d', ...
        animalID, sigPairs, upperPairs));
    set(gca, 'TickDir', 'out');
    box off;
end
end
