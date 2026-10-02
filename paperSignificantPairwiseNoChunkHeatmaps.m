function paperSignificantPairwiseNoChunkHeatmaps(combinedMatFile)
% Paper-formatted heatmaps of significant pairwise NO-CHUNK peak correlations

% Significance from  R.actualPairStats.sigFDRMask (on quest)

% The heatmap is reconstructed using:
%   actualPairStats.rows
%   actualPairStats.cols
%   actualPairStats.realVals

% Nonsignificant pairs remain NaN and are rendered black.

% This function DOES NOT:
%   - recompute significance using the old 3-SD shift-control criterion
%   - use Bayes factors
%   - use H0 / H50 / Hneg50

% D043 is excluded.

% RUN: paperSignificantPairwiseNoChunkXCorrHeatmaps

arguments
    combinedMatFile (1,1) string = ...
        "C:\Users\mirilab\Documents\GlobusTransfer\pairwiseNoChunkLagImbalanceBayes50ms_StoreyCorrThresh_DAVID_UPDATE.mat"
end

%% ========================================================================
%% LOAD RESULTS
%% ========================================================================

S = load(combinedMatFile, 'results');

if ~isfield(S, 'results')
    error('File does not contain variable "results":\n%s', combinedMatFile);
end

results = S.results;

if ~isfield(results, 'sessions') || isempty(results.sessions)
    error('results.sessions is missing or empty.');
end

sessions = results.sessions;

%% ========================================================================
%% EXPORT DIRECTORY
%% ========================================================================

saveDir = ...
    "X:\David\AnalysesData\InterneuronAnalyses\AA Paper Plots\pairwise no-chunk";

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% ========================================================================
%% FORMATTING
%% ========================================================================

titleFont = 10;
labelFont = 10;
tickFont = 10;
colorbarFont = 10;
axesLineWidth = 0.5;

%% ========================================================================
%% EXTRACT SESSIONS + REMOVE D043
%% ========================================================================

animalIDs = strings(0,1);
sessionData = cell(0,1);

for s = 1:numel(sessions)

    if iscell(sessions)
        R = sessions{s};
    else
        R = sessions(s);
    end

    if isempty(R)
        continue;
    end

    if isfield(R, 'animalID') && ~isempty(R.animalID)

        animalID = string(R.animalID);

    elseif isfield(R, 'baseDir') && ~isempty(R.baseDir)

        hit = regexp(char(R.baseDir), ...
            'D\d+', ...
            'match', ...
            'once');

        if isempty(hit)
            animalID = "Animal_" + s;
        else
            animalID = string(hit);
        end

    else
        animalID = "Animal_" + s;
    end

    % Exclude D043
    if strcmpi(animalID, "D043")
        continue;
    end

    animalIDs(end+1,1) = animalID; %#ok<AGROW>
    sessionData{end+1,1} = R; %#ok<AGROW>
end

nSess = numel(sessionData);

if nSess == 0
    error('No sessions remain after filtering.');
end

%% ========================================================================
%% CONSTRUCT SIGNIFICANT CORRELATION MATRICES
%% ========================================================================

sigMatCell = cell(nSess,1);

sigPairs = zeros(nSess,1);
totalPairs = zeros(nSess,1);

for s = 1:nSess

    R = sessionData{s};

    if ~isfield(R, 'actualPairStats') || ...
            isempty(R.actualPairStats)

        warning( ...
            '%s: actualPairStats missing or empty.', ...
            animalIDs(s));

        continue;
    end

    P = R.actualPairStats;

    requiredFields = { ...
        'rows', ...
        'cols', ...
        'realVals', ...
        'sigFDRMask'};

    for f = 1:numel(requiredFields)

        if ~isfield(P, requiredFields{f})

            error( ...
                '%s: actualPairStats.%s is missing.', ...
                animalIDs(s), ...
                requiredFields{f});
        end
    end

    rows = P.rows(:);
    cols = P.cols(:);
    realVals = P.realVals(:);
    sigMask = logical(P.sigFDRMask(:));

    if numel(rows) ~= numel(cols) || ...
            numel(rows) ~= numel(realVals) || ...
            numel(rows) ~= numel(sigMask)

        error( ...
            '%s: pairwise vectors have inconsistent lengths.', ...
            animalIDs(s));
    end

    %% ---- determine matrix dimensions ----

    if isfield(R, 'nInt') && ...
            ~isempty(R.nInt)

        nInt = double(R.nInt);

    else

        nInt = max(rows);
    end

    if isfield(R, 'nPyr') && ...
            ~isempty(R.nPyr)

        nPyr = double(R.nPyr);

    else

        nPyr = max(cols);
    end

    sigMat = nan(nInt, nPyr);

    %% ---- optional consistency check ----

    if isfield(R, 'actualSigMask') && ...
            ~isempty(R.actualSigMask)

        topMask = logical(R.actualSigMask(:));

        if numel(topMask) == numel(sigMask) && ...
                ~isequal(topMask, sigMask)

            warning( ...
                '%s: actualSigMask and actualPairStats.sigFDRMask differ.', ...
                animalIDs(s));
        end
    end

    %% ---- populate significant entries only ----

    for k = 1:numel(realVals)

        if sigMask(k) && ...
                isfinite(realVals(k)) && ...
                isfinite(rows(k)) && ...
                isfinite(cols(k))

            sigMat(rows(k), cols(k)) = realVals(k);
        end
    end

    sigMatCell{s} = sigMat;

    totalPairs(s) = numel(realVals);

    sigPairs(s) = sum( ...
        sigMask & ...
        isfinite(realVals));

    fprintf('\n=== %s ===\n', animalIDs(s));
    fprintf('Pair matrix: %d interneurons x %d pyramidal neurons\n', ...
        nInt, nPyr);
    fprintf('Tested pairs: %d\n', totalPairs(s));
    fprintf('Storey/FDR-significant pairs: %d\n', sigPairs(s));

    if isfield(P, 'qAlpha')
        fprintf('qAlpha: %.6g\n', P.qAlpha);
    end

    if isfield(P, 'corrThresh')
        fprintf('corrThresh: %.6g\n', P.corrThresh);
    end
end

%% ========================================================================
%% COMMON COLOR LIMITS
%% ========================================================================

allSigCorrs = [];

for s = 1:nSess

    sigMat = sigMatCell{s};

    if isempty(sigMat)
        continue;
    end

    vals = sigMat(isfinite(sigMat));

    allSigCorrs = [ ...
        allSigCorrs; ...
        vals(:)]; %#ok<AGROW>
end

if isempty(allSigCorrs)

    commonCLim = [0 1];

else

    cMin = min(allSigCorrs);
    cMax = max(allSigCorrs);

    if cMin == cMax
        cPad = 0.01;
    else
        cPad = 0.05 * (cMax - cMin);
    end

    commonCLim = [ ...
        cMin - cPad, ...
        cMax + cPad];
end

%% ========================================================================
%% TILED HEATMAP FIGURE
%% ========================================================================

fig = figure( ...
    'Name', 'Pairwise no-chunk significant peak correlations', ...
    'Color', 'w', ...
    'Position', [100 100 1500 650]);

tileLay = tiledlayout( ...
    1, nSess, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

title(tileLay, ...
    'Pairwise no-chunk significant peak correlations', ...
    'FontSize', titleFont);

for s = 1:nSess

    ax = nexttile(tileLay, s);

    sigMat = sigMatCell{s};

    if isempty(sigMat)

        title(ax, ...
            sprintf('%s (missing)', animalIDs(s)), ...
            'FontSize', titleFont);

        axis(ax, 'off');
        continue;
    end

    hImg = imagesc(ax, sigMat);

    % Significant entries are opaque.
    % Nonsignificant NaNs are transparent.
    set(hImg, ...
        'AlphaData', ...
        ~isnan(sigMat));

    % Black background shows through nonsignificant entries.
    set(ax, ...
        'Color', 'k');

    axis(ax, 'xy');

    colormap(ax, parula);
    clim(ax, commonCLim);

    xlabel(ax, ...
        'Pyramidal neurons', ...
        'FontSize', labelFont);

    ylabel(ax, ...
        'Interneurons', ...
        'FontSize', labelFont);

    title(ax, ...
        sprintf('%s (%d / %d)', ...
        animalIDs(s), ...
        sigPairs(s), ...
        totalPairs(s)), ...
        'FontSize', titleFont);

    box(ax, 'off');

    set(ax, ...
        'FontSize', tickFont, ...
        'LineWidth', axesLineWidth, ...
        'TickDir', 'out', ...
        'TickLength', [0.025 0.025], ...
        'XColor', 'k', ...
        'YColor', 'k');
end

%% ========================================================================
%% SHARED COLORBAR
%% ========================================================================

c = colorbar;
c.Layout.Tile = 'east';

ylabel(c, ...
    'Peak correlation', ...
    'FontSize', labelFont);

c.FontSize = colorbarFont;

%% ========================================================================
%% EXPORT
%% ========================================================================

exportgraphics(fig, ...
    fullfile(saveDir, ...
    'pairwise_nochunk_significant_peak_corr_heatmaps.pdf'), ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');

fprintf('\nPairwise no-chunk heatmap saved to:\n%s\n', saveDir);

end
