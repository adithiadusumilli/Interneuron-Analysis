function plotAllVsAllVsIntPyrBias_FDR(fdrMatFile)
% compares all-vs-all pair distributions to int-vs-pyr distributions
% to see whether int-vs-pyr pairs show a bias relative to the full pair population

% works for either NO-CHUNK or CHUNKED FDR files

S = load(fdrMatFile, 'FDRresults');
FDRresults = S.FDRresults;

numSessions = numel(FDRresults.sessions);

for sess = 1:numSessions
    D = FDRresults.sessions{sess};

    % detect file type
    if isfield(D, 'peakCorrMatAll')
        peakCorrs = D.peakCorrMatAll;
        peakLags = D.peakLagMatAll;

        if isfield(D, 'typeVec') && ~isempty(D.typeVec)
            typeVec = D.typeVec(:);
            intIdx = find(typeVec == 1);
            pyrIdx = find(typeVec == 0);
        else
            warning('sess %d: no typeVec found, skipping', sess);
            continue;
        end

    elseif isfield(D, 'peakCorrMat_all')
        peakCorrs = D.peakCorrMat_all;
        peakLags = D.peakLagSecMat_all;

        if isfield(D, 'intIdx') && ~isempty(D.intIdx)
            intIdx = D.intIdx(:);
        elseif isfield(D, 'neuronType') && ~isempty(D.neuronType)
            intIdx = find(D.neuronType == 1);
        else
            warning('sess %d: no intIdx found, skipping', sess);
            continue;
        end

        if isfield(D, 'pyrIdx') && ~isempty(D.pyrIdx)
            pyrIdx = D.pyrIdx(:);
        elseif isfield(D, 'neuronType') && ~isempty(D.neuronType)
            pyrIdx = find(D.neuronType == 0);
        else
            warning('sess %d: no pyrIdx found, skipping', sess);
            continue;
        end
    else
        warning('sess %d: unrecognized file structure, skipping', sess);
        continue;
    end

    if isempty(intIdx) || isempty(pyrIdx)
        warning('sess %d: empty intIdx or pyrIdx, skipping', sess);
        continue;
    end

    animalID = D.animalID;
    if isempty(animalID)
        animalID = sprintf('sess %d', sess);
    end

    % ---- all-vs-all upper triangle ----
    upperMask = triu(true(size(peakCorrs)), 1);
    allCorrVec = peakCorrs(upperMask);
    allLagVec = peakLags(upperMask);

    % ---- int-vs-pyr block ----
    intPyrCorr = peakCorrs(intIdx, pyrIdx);
    intPyrLag = peakLags(intIdx, pyrIdx);

    intPyrCorrVec = intPyrCorr(:);
    intPyrLagVec = intPyrLag(:);

    % ---- lag bins ----
    lagAll = [allLagVec(:); intPyrLagVec(:)];
    lagAll = lagAll(~isnan(lagAll));

    if isempty(lagAll)
        minLagSec = -0.2;
        maxLagSec = 0.2;
    else
        minLagSec = floor(min(lagAll)*1000)/1000;
        maxLagSec = ceil(max(lagAll)*1000)/1000;
    end
    lagBinEdgesSec = (minLagSec-0.001):0.001:(maxLagSec+0.001);

    % ==================== peak correlation comparison ====================
    figure('Name', sprintf('%s – all-vs-all vs int-pyr peak corr', animalID), 'Color', 'w');
    hold on;
    histogram(allCorrVec, 50, 'Normalization', 'probability', 'FaceAlpha', 0.45, 'EdgeColor', 'none');
    histogram(intPyrCorrVec, 50, 'Normalization', 'probability', 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    xlabel('peak correlation');
    ylabel('probability');
    title(sprintf('%s – all-vs-all vs int-pyr peak correlations', animalID));
    legend({'all-vs-all','int-vs-pyr'}, 'Location', 'best');
    grid on;

    % ==================== peak lag comparison ====================
    figure('Name', sprintf('%s – all-vs-all vs int-pyr peak lag', animalID), 'Color', 'w');
    hold on;
    histogram(allLagVec, 'BinEdges', lagBinEdgesSec, 'Normalization', 'probability', 'FaceAlpha', 0.45, 'EdgeColor', 'none');
    histogram(intPyrLagVec, 'BinEdges', lagBinEdgesSec, 'Normalization', 'probability', 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    xlabel('peak lag (s)');
    ylabel('probability');
    title(sprintf('%s – all-vs-all vs int-pyr peak lags', animalID));
    legend({'all-vs-all','int-vs-pyr'}, 'Location', 'best');
    grid on;

    % ==================== scatter comparison ====================
    figure('Name', sprintf('%s – all-vs-all vs int-pyr lag vs corr', animalID), 'Color', 'w');
    hold on;
    scatter(allLagVec(:), allCorrVec(:), 8, '.', 'MarkerEdgeAlpha', 0.25);
    scatter(intPyrLagVec(:), intPyrCorrVec(:), 12, '.');
    xlabel('peak lag (s)');
    ylabel('peak correlation');
    title(sprintf('%s – all-vs-all vs int-pyr lag vs corr', animalID));
    legend({'all-vs-all','int-vs-pyr'}, 'Location', 'best');
    grid on;

    fprintf('\n%s\n', animalID);
    fprintf('all-vs-all pairs: %d\n', numel(allCorrVec));
    fprintf('int-vs-pyr pairs: %d\n', numel(intPyrCorrVec));
end
end
