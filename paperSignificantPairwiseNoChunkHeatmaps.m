function paperSignificantPairwiseNoChunkHeatmaps(combinedMatFile)
% Heatmaps of significant pairwise NO-CHUNK peak correlations per session.
%
% Uses new Storey/FDR-corrected pairwise output:
%   R.actualPairStats.rows
%   R.actualPairStats.cols
%   R.actualPairStats.realVals
%   R.actualPairStats.sigFDRMask
%
% rows = local interneuron indices (1:nInt)
% cols = global neuron indices (nInt+1:nAll)
%
% Therefore, pyramidal-neuron column indices are converted to:
%   localPyrCol = cols - nInt
%
% Significance is taken directly from the saved Storey/FDR result.
% No 3-SD shift-control significance is recomputed.
% No Bayes-factor variables are used.
%
% Nonsignificant entries are set to NaN and rendered transparent.
%
% D043 is excluded.
%
% RUN:
%   paperSignificantPairwiseNoChunkHeatmaps
%
% Or:
%   paperSignificantPairwiseNoChunkHeatmaps(combinedMatFile)

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

%% ---- loop through sessions ----

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

    if ~isfield(P, 'rows') || ...
       ~isfield(P, 'cols') || ...
       ~isfield(P, 'realVals') || ...
       ~isfield(P, 'sigFDRMask')

        warning( ...
            'sess %d: required actualPairStats fields missing, skipping', ...
            sess);

        continue;
    end

    rows = P.rows(:);
    cols = P.cols(:);
    realVals = P.realVals(:);

    % Storey/FDR-corrected significance calculated upstream
    sigMask = logical(P.sigFDRMask(:));

    if numel(rows) ~= numel(cols) || ...
       numel(rows) ~= numel(realVals) || ...
       numel(rows) ~= numel(sigMask)

        warning( ...
            'sess %d: pairwise vectors have inconsistent lengths, skipping', ...
            sess);

        continue;
    end

    %% ---- matrix dimensions ----

    nInt = double(R.nInt);
    nPyr = double(R.nPyr);

    %% ---- convert global pyramidal indices to local pyramidal indices ----
    %
    % Example D020:
    %
    %   interneurons = global neurons 1:15
    %   pyramidal    = global neurons 16:82
    %
    % For the nInt x nPyr heatmap:
    %
    %   global pyr 16 -> local pyr 1
    %   global pyr 17 -> local pyr 2
    %   ...
    %   global pyr 82 -> local pyr 67

    pyrCols = cols - nInt;

    %% ---- sanity check indexing ----

    if any(rows < 1 | rows > nInt)

        error( ...
            'sess %d (%s): interneuron row indices fall outside 1:nInt.', ...
            sess, animalID);
    end

    if any(pyrCols < 1 | pyrCols > nPyr)

        error( ...
            'sess %d (%s): converted pyramidal indices fall outside 1:nPyr.', ...
            sess, animalID);
    end

    %% ---- optional check against top-level significance mask ----

    if isfield(R, 'actualSigMask') && ~isempty(R.actualSigMask)

        actualSigMask = logical(R.actualSigMask(:));

        if numel(actualSigMask) == numel(sigMask) && ...
                ~isequal(actualSigMask, sigMask)

            warning( ...
                'sess %d (%s): actualSigMask differs from sigFDRMask', ...
                sess, animalID);
        end
    end

    %% ---- construct significant correlation matrix ----

    sigMat = nan(nInt, nPyr);
    sigPairs = 0;

    for k = 1:numel(realVals)

        if sigMask(k) && ...
           ~isnan(realVals(k)) && ...
           ~isnan(rows(k)) && ...
           ~isnan(pyrCols(k))

            sigMat(rows(k), pyrCols(k)) = realVals(k);

            sigPairs = sigPairs + 1;
        end
    end

    totalPairs = nInt * nPyr;

    %% ---- print information ----

    fprintf('\n=== sess %d: %s ===\n', sess, animalID);
    fprintf('matrix size: %d interneurons x %d pyramidal neurons\n', ...
        nInt, nPyr);
    fprintf('pairs: %d total | %d significant\n', ...
        totalPairs, sigPairs);

    %% ---- base label ----

    baseLabel = sprintf('sess %d', sess);

    if isfield(R, 'baseDir') && ~isempty(R.baseDir)

        baseLabel = sprintf( ...
            'sess %d – %s', ...
            sess, ...
            R.baseDir);
    end

    %% ---- plot heatmap ----

    figure( ...
        'Name', ...
        sprintf('sess %d – NO-CHUNK significant peak corr', sess), ...
        'Color', 'w');

    hImg = imagesc(sigMat);

    set(hImg, ...
        'AlphaData', ...
        ~isnan(sigMat));  % hide nonsig pairs

    set(gca, ...
        'Color', 'k');    % black background for NaNs

    axis xy;

    colormap(parula);

    c = colorbar;
    ylabel(c, 'peak correlation');

    xlabel('pyramidal neurons');
    ylabel('interneurons');

    title(sprintf( ...
        '%s – NO-CHUNK significant pairs: %d / %d', ...
        baseLabel, ...
        sigPairs, ...
        totalPairs));

    set(gca, 'TickDir', 'out');

    box off;

end

end
