function paperSignificantPairwiseNoChunkHeatmaps(combinedMatFile)
% Heatmaps of significant pairwise NO-CHUNK peak correlations per session.
%
% Uses new Storey/FDR-corrected pairwise output:
%   R.actualPairStats.rows
%   R.actualPairStats.cols
%   R.actualPairStats.realVals
%   R.actualPairStats.sigFDRMask
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

    if isfield(R, 'nInt') && ~isempty(R.nInt)
        nInt = double(R.nInt);
    else
        nInt = max(rows);
    end

    if isfield(R, 'nPyr') && ~isempty(R.nPyr)
        nPyr = double(R.nPyr);
    else
        nPyr = max(cols);
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
           ~isnan(cols(k))

            sigMat(rows(k), cols(k)) = realVals(k);
            sigPairs = sigPairs + 1;
        end
    end

    % Same meaning as old nInt*nPyr matrix.
    totalPairs = nInt * nPyr;

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

    % Hide nonsignificant pairs
    set(hImg, ...
        'AlphaData', ...
        ~isnan(sigMat));

    % Black background for nonsignificant/NaN entries
    set(gca, 'Color', 'k');

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
