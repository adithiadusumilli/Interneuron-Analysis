function paperSignificantPairwiseNoChunk(combinedMatFile)
% Plot significant pairwise NO-CHUNK xcorrs per animal using new
% Storey/FDR-corrected pairwise output.
%
% Uses:
%   R.actualPairStats.realVals      pairwise peak correlations
%   R.actualPairStats.lagVals       pairwise peak lags (seconds)
%   R.actualPairStats.sigFDRMask    Storey/FDR significance mask
%
% Significance is taken directly from the saved Storey/FDR result.
% No 3-SD shift-control significance is recomputed.
% No Bayes-factor variables are used.
%
% Also plots:
%   - peak corr > 0.2 (no significance)
%   - significant |lag| > 0.2 s
%   - scatterhist of significant pairs
%   - peak corr >= 0.2 (no significance)
%   - peak lags corresponding to corr >= 0.2
%   - scatter of lag vs corr for corr >= 0.2
%
% D043 is excluded.
%
% RUN:
%   paperSignificantPairwiseNoChunk
%
% Or:
%   paperSignificantPairwiseNoChunk(combinedMatFile)

arguments
    combinedMatFile (1,1) string = ...
        "C:\Users\mirilab\Documents\GlobusTransfer\pairwiseNoChunkLagImbalanceBayes50ms_StoreyCorrThresh_DAVID_UPDATE.mat"
end

%% ---- load results ----

S = load(combinedMatFile, 'results');

if ~isfield(S, 'results')
    error('File does not contain variable "results":\n%s', combinedMatFile);
end

results = S.results;

if ~isfield(results, 'sessions') || isempty(results.sessions)
    error('results.sessions is missing or empty.');
end

sessions = results.sessions;
numSessions = numel(sessions);

%% ---- loop through animals ----

for sess = 1:numSessions

    if iscell(sessions)
        R = sessions{sess};
    else
        R = sessions(sess);
    end

    if isempty(R)
        fprintf('sess %d: empty data, skipping\n', sess);
        continue;
    end

    %% ---- animal ID ----

    if isfield(R, 'animalID') && ~isempty(R.animalID)

        animalID = string(R.animalID);

    elseif isfield(R, 'baseDir') && ~isempty(R.baseDir)

        hit = regexp(char(R.baseDir), 'D\d+', 'match', 'once');

        if isempty(hit)
            animalID = "Animal_" + sess;
        else
            animalID = string(hit);
        end

    else
        animalID = "Animal_" + sess;
    end

    %% ---- exclude D043 ----

    if strcmpi(animalID, "D043")
        fprintf('sess %d (%s): excluded\n', sess, animalID);
        continue;
    end

    %% ---- get pairwise data ----

    if ~isfield(R, 'actualPairStats') || isempty(R.actualPairStats)
        fprintf('sess %d: actualPairStats empty, skipping\n', sess);
        continue;
    end

    P = R.actualPairStats;

    if ~isfield(P, 'realVals') || ...
       ~isfield(P, 'lagVals') || ...
       ~isfield(P, 'sigFDRMask')

        warning( ...
            'sess %d: required actualPairStats fields missing, skipping', ...
            sess);

        continue;
    end

    allCorrVec = P.realVals(:);
    allLagVec = P.lagVals(:);

    % Storey/FDR-corrected significance calculated upstream
    sigMask = logical(P.sigFDRMask(:));

    if numel(allCorrVec) ~= numel(allLagVec) || ...
            numel(allCorrVec) ~= numel(sigMask)

        warning( ...
            'sess %d: size mismatch between pairwise vectors, skipping', ...
            sess);

        continue;
    end

    %% ---- optional check against saved top-level significance mask ----

    if isfield(R, 'actualSigMask') && ~isempty(R.actualSigMask)

        actualSigMask = logical(R.actualSigMask(:));

        if numel(actualSigMask) == numel(sigMask) && ...
                ~isequal(actualSigMask, sigMask)

            warning( ...
                'sess %d (%s): actualSigMask differs from sigFDRMask', ...
                sess, animalID);
        end
    end

    %% ---- significant pairs ----

    sigCorrVec = allCorrVec(sigMask);
    sigLagVec = allLagVec(sigMask);

    totalPairs = numel(allCorrVec);
    sigPairs = numel(sigCorrVec);

    fprintf('\n=== sess %d: %s ===\n', sess, animalID);

    if isfield(R, 'baseDir') && ~isempty(R.baseDir)
        fprintf('baseDir: %s\n', R.baseDir);
    end

    fprintf( ...
        'pairs: %d total | %d significant\n', ...
        totalPairs, sigPairs);

    if isfield(P, 'qAlpha') && ~isempty(P.qAlpha)
        fprintf('Storey/FDR qAlpha: %.6g\n', P.qAlpha);
    end

    if isfield(P, 'corrThresh') && ~isempty(P.corrThresh)
        fprintf('correlation threshold: %.6g\n', P.corrThresh);
    end

    %% ---- 1 ms bin edges for lag histogram (sec) ----

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

    lagBinEdgesSec = ...
        (minLagSec-0.001):0.001:(maxLagSec+0.001);

    %% ---- plot: significant peak lags ----

    figure;

    histogram( ...
        sigLagVec, ...
        'BinEdges', lagBinEdgesSec, ...
        'FaceAlpha', 0.8, ...
        'EdgeColor', 'none');

    xlabel('peak lag (s)');
    ylabel('count');

    title(sprintf( ...
        'sess %d – NO-CHUNK significant peak lags (n=%d / %d)', ...
        sess, ...
        sigPairs, ...
        totalPairs));

    grid on;

    %% ---- plot: significant peak correlations ----

    figure;

    histogram( ...
        sigCorrVec, ...
        50, ...
        'FaceAlpha', 0.8, ...
        'EdgeColor', 'none');

    xlabel('peak correlation');
    ylabel('count');

    title(sprintf( ...
        'sess %d – NO-CHUNK significant peak correlations (n=%d / %d)', ...
        sess, ...
        sigPairs, ...
        totalPairs));

    grid on;

    %% ---- plot: peak correlations > 0.2 (no significance) ----

    corrOverThresh = allCorrVec(allCorrVec > 0.2);

    figure;

    histogram( ...
        corrOverThresh, ...
        50, ...
        'FaceAlpha', 0.8, ...
        'EdgeColor', 'none');

    xlabel('peak correlation (> 0.2)');
    ylabel('count');

    title(sprintf( ...
        'sess %d – NO-CHUNK peak correlations > 0.2 (n=%d / %d)', ...
        sess, ...
        numel(corrOverThresh), ...
        totalPairs));

    grid on;

    %% ---- plot: significant |lag| > 0.2 s ----

    lagsOverThresh = ...
        sigLagVec(abs(sigLagVec) > 0.2);

    figure;

    histogram( ...
        lagsOverThresh, ...
        'BinEdges', lagBinEdgesSec, ...
        'FaceAlpha', 0.8, ...
        'EdgeColor', 'none');

    xlabel('peak lag (s)');
    ylabel('count');

    title(sprintf( ...
        'sess %d – NO-CHUNK significant |peak lag| > 0.2 s (n=%d / %d sig)', ...
        sess, ...
        numel(lagsOverThresh), ...
        numel(sigLagVec)));

    grid on;

    %% ---- scatterhist: significant pairs ----

    if ~isempty(sigLagVec)

        figure;

        scatterhist( ...
            sigLagVec(:), ...
            sigCorrVec(:), ...
            'Direction', 'out', ...
            'Marker', '.');

        xlabel('peak lag (s)');
        ylabel('peak correlation');

        title(sprintf( ...
            'sess %d – NO-CHUNK significant pairs (n=%d / %d)', ...
            sess, ...
            sigPairs, ...
            totalPairs));
    end

    %% ---- plot: peak correlations >= 0.2 (no significance) ----

    corrMask = ...
        (allCorrVec >= 0.2) & ...
        ~isnan(allCorrVec) & ...
        ~isnan(allLagVec);

    corrOverThresh = allCorrVec(corrMask);
    lagsForCorrOverThresh = allLagVec(corrMask);

    figure;

    histogram( ...
        corrOverThresh, ...
        50, ...
        'FaceAlpha', 0.8, ...
        'EdgeColor', 'none');

    xlabel('peak correlation (>= 0.2)');
    ylabel('count');

    title(sprintf( ...
        'sess %d – NO-CHUNK peak correlations >= 0.2 (n=%d / %d)', ...
        sess, ...
        numel(corrOverThresh), ...
        totalPairs));

    grid on;

    %% ---- peak LAGS corresponding to corr >= 0.2 ----

    if ~isempty(lagsForCorrOverThresh)

        % 1 ms bin edges (sec) for these lags
        minLagSec2 = ...
            floor(min(lagsForCorrOverThresh)*1000)/1000;

        maxLagSec2 = ...
            ceil(max(lagsForCorrOverThresh)*1000)/1000;

        lagBinEdgesSec2 = ...
            (minLagSec2-0.001):0.001:(maxLagSec2+0.001);

        figure;

        histogram( ...
            lagsForCorrOverThresh, ...
            'BinEdges', lagBinEdgesSec2, ...
            'FaceAlpha', 0.8, ...
            'EdgeColor', 'none');

        xlabel('peak lag (s) for pairs with corr >= 0.2');
        ylabel('count');

        title(sprintf( ...
            'sess %d – NO-CHUNK peak lags where peak corr >= 0.2 (n=%d)', ...
            sess, ...
            numel(lagsForCorrOverThresh)));

        grid on;

        %% ---- scatter lag vs corr for corr >= 0.2 ----

        figure;

        scatter( ...
            lagsForCorrOverThresh(:), ...
            corrOverThresh(:), ...
            '.');

        xlabel('peak lag (s)');
        ylabel('peak correlation');

        title(sprintf( ...
            'sess %d – NO-CHUNK lag vs corr for corr >= 0.2 (n=%d)', ...
            sess, ...
            numel(corrOverThresh)));

        grid on;
    end
end

end
