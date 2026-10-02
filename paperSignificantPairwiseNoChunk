function paperSignificantPairwiseNoChunk(combinedMatFile)
% Paper-formatted pairwise NO-CHUNK cross-correlation plots.

% Significance is taken DIRECTLY from the saved Storey/FDR-corrected
% significance mask:

%   R.actualPairStats.sigFDRMask

% This function DOES NOT:
%   - recompute significance using the old 3-SD shift-control criterion
%   - use Bayes factors
%   - use H0 / H50 / Hneg50

% Uses:
%   actualPairStats.rows
%   actualPairStats.cols
%   actualPairStats.realVals
%   actualPairStats.lagVals
%   actualPairStats.sigFDRMask

% D043 is excluded!

% RUN: paperSignificantPairwiseNoChunkXCorrsPerAnimal

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
axesLineWidth = 0.5;

histColor = [0.3 0.6 0.8];
scatterColor = [0 0 0];

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

        hit = regexp(char(R.baseDir), 'D\d+', 'match', 'once');

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
%% EXTRACT PAIRWISE DATA
%% ========================================================================

allCorrCell = cell(nSess,1);
allLagCell = cell(nSess,1);

sigCorrCell = cell(nSess,1);
sigLagCell = cell(nSess,1);

corrOverCell = cell(nSess,1);
lagCorrOverCell = cell(nSess,1);

sigLagOverCell = cell(nSess,1);

totalPairs = zeros(nSess,1);
sigPairs = zeros(nSess,1);

for s = 1:nSess

    R = sessionData{s};

    if ~isfield(R, 'actualPairStats') || isempty(R.actualPairStats)
        warning('%s: actualPairStats missing or empty.', animalIDs(s));
        continue;
    end

    P = R.actualPairStats;

    requiredFields = { ...
        'realVals', ...
        'lagVals', ...
        'sigFDRMask'};

    for f = 1:numel(requiredFields)
        if ~isfield(P, requiredFields{f})
            error( ...
                '%s: actualPairStats.%s is missing.', ...
                animalIDs(s), ...
                requiredFields{f});
        end
    end

    allCorrVec = P.realVals(:);
    allLagVec = P.lagVals(:);
    sigMask = logical(P.sigFDRMask(:));

    if numel(allCorrVec) ~= numel(allLagVec) || ...
            numel(allCorrVec) ~= numel(sigMask)

        error( ...
            '%s: realVals, lagVals, and sigFDRMask have different lengths.', ...
            animalIDs(s));
    end

    % Optional consistency check against saved top-level mask
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

    validMask = ...
        ~isnan(allCorrVec) & ...
        ~isnan(allLagVec);

    sigValidMask = ...
        sigMask & validMask;

    sigCorrVec = allCorrVec(sigValidMask);
    sigLagVec = allLagVec(sigValidMask);

    % Correlation threshold plots are NOT restricted to significant pairs.
    corrMask = ...
        validMask & ...
        allCorrVec >= 0.2;

    corrOverVec = allCorrVec(corrMask);
    lagCorrOverVec = allLagVec(corrMask);

    sigLagOverVec = ...
        sigLagVec(abs(sigLagVec) > 0.2);

    allCorrCell{s} = allCorrVec(validMask);
    allLagCell{s} = allLagVec(validMask);

    sigCorrCell{s} = sigCorrVec;
    sigLagCell{s} = sigLagVec;

    corrOverCell{s} = corrOverVec;
    lagCorrOverCell{s} = lagCorrOverVec;

    sigLagOverCell{s} = sigLagOverVec;

    totalPairs(s) = sum(validMask);
    sigPairs(s) = sum(sigValidMask);

    fprintf('\n=== %s ===\n', animalIDs(s));
    fprintf('Valid pairs: %d\n', totalPairs(s));
    fprintf('Storey/FDR-significant pairs: %d\n', sigPairs(s));

    if isfield(P, 'qAlpha')
        fprintf('qAlpha: %.6g\n', P.qAlpha);
    end

    if isfield(P, 'corrThresh')
        fprintf('corrThresh: %.6g\n', P.corrThresh);
    end
end

%% ========================================================================
%% COMMON LIMITS: SIGNIFICANT PEAK LAGS
%% ========================================================================

allSigLags = vertcat(sigLagCell{:});
allSigLags = allSigLags(isfinite(allSigLags));

if isempty(allSigLags)

    commonSigLagXLim = [-0.2 0.2];

else

    xMin = floor(min(allSigLags) * 1000) / 1000;
    xMax = ceil(max(allSigLags) * 1000) / 1000;

    if xMin == xMax
        xMin = xMin - 0.001;
        xMax = xMax + 0.001;
    else
        xMin = xMin - 0.001;
        xMax = xMax + 0.001;
    end

    commonSigLagXLim = [xMin xMax];
end

sigLagEdges = ...
    commonSigLagXLim(1):0.001:commonSigLagXLim(2);

if numel(sigLagEdges) < 2
    sigLagEdges = linspace( ...
        commonSigLagXLim(1), ...
        commonSigLagXLim(2), ...
        3);
end

maxCount = 0;

for s = 1:nSess

    vals = sigLagCell{s};

    if isempty(vals)
        continue;
    end

    counts = histcounts(vals, sigLagEdges);

    if ~isempty(counts)
        maxCount = max(maxCount, max(counts));
    end
end

if maxCount == 0
    commonSigLagYLim = [0 1];
else
    commonSigLagYLim = [0 1.05 * maxCount];
end

%% ========================================================================
%% FIGURE 1: SIGNIFICANT PEAK LAGS
%% ========================================================================

fig1 = figure( ...
    'Name', 'Pairwise no-chunk significant peak lags', ...
    'Color', 'w', ...
    'Position', [100 100 1500 650]);

tile1 = tiledlayout( ...
    1, nSess, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

title(tile1, ...
    'Pairwise no-chunk significant peak lags', ...
    'FontSize', titleFont);

for s = 1:nSess

    ax = nexttile(tile1, s);
    hold(ax, 'on');

    vals = sigLagCell{s};

    if ~isempty(vals)

        histogram(ax, ...
            vals, ...
            'BinEdges', sigLagEdges, ...
            'FaceColor', histColor, ...
            'EdgeColor', 'none');
    end

    xlim(ax, commonSigLagXLim);
    ylim(ax, commonSigLagYLim);

    xlabel(ax, ...
        'Peak lag (s)', ...
        'FontSize', labelFont);

    ylabel(ax, ...
        'Count', ...
        'FontSize', labelFont);

    title(ax, ...
        sprintf('%s (n=%d / %d)', ...
        animalIDs(s), ...
        sigPairs(s), ...
        totalPairs(s)), ...
        'FontSize', titleFont);

    box(ax, 'off');
    grid(ax, 'off');

    set(ax, ...
        'FontSize', tickFont, ...
        'LineWidth', axesLineWidth, ...
        'TickDir', 'out', ...
        'TickLength', [0.025 0.025], ...
        'XColor', 'k', ...
        'YColor', 'k');
end

exportgraphics(fig1, ...
    fullfile(saveDir, ...
    'pairwise_nochunk_significant_peak_lags.pdf'), ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');

%% ========================================================================
%% COMMON LIMITS: SIGNIFICANT PEAK CORRELATIONS
%% ========================================================================

allSigCorrs = vertcat(sigCorrCell{:});
allSigCorrs = allSigCorrs(isfinite(allSigCorrs));

if isempty(allSigCorrs)

    commonSigCorrXLim = [0 1];

else

    xMin = min(allSigCorrs);
    xMax = max(allSigCorrs);

    if xMin == xMax
        xPad = 0.01;
    else
        xPad = 0.05 * (xMax - xMin);
    end

    commonSigCorrXLim = ...
        [xMin - xPad, xMax + xPad];
end

sigCorrEdges = linspace( ...
    commonSigCorrXLim(1), ...
    commonSigCorrXLim(2), ...
    51);

maxCount = 0;

for s = 1:nSess

    vals = sigCorrCell{s};

    if isempty(vals)
        continue;
    end

    counts = histcounts(vals, sigCorrEdges);

    if ~isempty(counts)
        maxCount = max(maxCount, max(counts));
    end
end

if maxCount == 0
    commonSigCorrYLim = [0 1];
else
    commonSigCorrYLim = [0 1.05 * maxCount];
end

%% ========================================================================
%% FIGURE 2: SIGNIFICANT PEAK CORRELATIONS
%% ========================================================================

fig2 = figure( ...
    'Name', 'Pairwise no-chunk significant peak correlations', ...
    'Color', 'w', ...
    'Position', [100 100 1500 650]);

tile2 = tiledlayout( ...
    1, nSess, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

title(tile2, ...
    'Pairwise no-chunk significant peak correlations', ...
    'FontSize', titleFont);

for s = 1:nSess

    ax = nexttile(tile2, s);
    hold(ax, 'on');

    vals = sigCorrCell{s};

    if ~isempty(vals)

        histogram(ax, ...
            vals, ...
            'BinEdges', sigCorrEdges, ...
            'FaceColor', histColor, ...
            'EdgeColor', 'none');
    end

    xlim(ax, commonSigCorrXLim);
    ylim(ax, commonSigCorrYLim);

    xlabel(ax, ...
        'Peak correlation', ...
        'FontSize', labelFont);

    ylabel(ax, ...
        'Count', ...
        'FontSize', labelFont);

    title(ax, ...
        sprintf('%s (n=%d / %d)', ...
        animalIDs(s), ...
        sigPairs(s), ...
        totalPairs(s)), ...
        'FontSize', titleFont);

    box(ax, 'off');
    grid(ax, 'off');

    set(ax, ...
        'FontSize', tickFont, ...
        'LineWidth', axesLineWidth, ...
        'TickDir', 'out', ...
        'TickLength', [0.025 0.025], ...
        'XColor', 'k', ...
        'YColor', 'k');
end

exportgraphics(fig2, ...
    fullfile(saveDir, ...
    'pairwise_nochunk_significant_peak_correlations.pdf'), ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');

%% ========================================================================
%% FIGURE 3: SIGNIFICANT |PEAK LAG| > 0.2 S
%% ========================================================================

maxCount = 0;

for s = 1:nSess

    vals = sigLagOverCell{s};

    if isempty(vals)
        continue;
    end

    counts = histcounts(vals, sigLagEdges);

    if ~isempty(counts)
        maxCount = max(maxCount, max(counts));
    end
end

if maxCount == 0
    commonLagOverYLim = [0 1];
else
    commonLagOverYLim = [0 1.05 * maxCount];
end

fig3 = figure( ...
    'Name', 'Pairwise no-chunk significant absolute peak lags over 0.2 s', ...
    'Color', 'w', ...
    'Position', [100 100 1500 650]);

tile3 = tiledlayout( ...
    1, nSess, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

title(tile3, ...
    'Pairwise no-chunk significant |peak lag| > 0.2 s', ...
    'FontSize', titleFont);

for s = 1:nSess

    ax = nexttile(tile3, s);
    hold(ax, 'on');

    vals = sigLagOverCell{s};

    if ~isempty(vals)

        histogram(ax, ...
            vals, ...
            'BinEdges', sigLagEdges, ...
            'FaceColor', histColor, ...
            'EdgeColor', 'none');
    end

    xlim(ax, commonSigLagXLim);
    ylim(ax, commonLagOverYLim);

    xlabel(ax, ...
        'Peak lag (s)', ...
        'FontSize', labelFont);

    ylabel(ax, ...
        'Count', ...
        'FontSize', labelFont);

    title(ax, ...
        sprintf('%s (n=%d / %d significant)', ...
        animalIDs(s), ...
        numel(vals), ...
        sigPairs(s)), ...
        'FontSize', titleFont);

    box(ax, 'off');
    grid(ax, 'off');

    set(ax, ...
        'FontSize', tickFont, ...
        'LineWidth', axesLineWidth, ...
        'TickDir', 'out', ...
        'TickLength', [0.025 0.025], ...
        'XColor', 'k', ...
        'YColor', 'k');
end

exportgraphics(fig3, ...
    fullfile(saveDir, ...
    'pairwise_nochunk_significant_peak_lags_over_0p2s.pdf'), ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');

%% ========================================================================
%% FIGURE 4: SIGNIFICANT LAG VS CORRELATION
%% ========================================================================

fig4 = figure( ...
    'Name', 'Pairwise no-chunk significant lag vs correlation', ...
    'Color', 'w', ...
    'Position', [100 100 1500 650]);

tile4 = tiledlayout( ...
    1, nSess, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

title(tile4, ...
    'Pairwise no-chunk significant lag vs. correlation', ...
    'FontSize', titleFont);

for s = 1:nSess

    ax = nexttile(tile4, s);
    hold(ax, 'on');

    lagVals = sigLagCell{s};
    corrVals = sigCorrCell{s};

    if ~isempty(lagVals)

        scatter(ax, ...
            lagVals, ...
            corrVals, ...
            12, ...
            scatterColor, ...
            'filled');
    end

    xlim(ax, commonSigLagXLim);
    ylim(ax, commonSigCorrXLim);

    xlabel(ax, ...
        'Peak lag (s)', ...
        'FontSize', labelFont);

    ylabel(ax, ...
        'Peak correlation', ...
        'FontSize', labelFont);

    title(ax, ...
        sprintf('%s (n=%d)', ...
        animalIDs(s), ...
        sigPairs(s)), ...
        'FontSize', titleFont);

    box(ax, 'off');
    grid(ax, 'off');

    set(ax, ...
        'FontSize', tickFont, ...
        'LineWidth', axesLineWidth, ...
        'TickDir', 'out', ...
        'TickLength', [0.025 0.025], ...
        'XColor', 'k', ...
        'YColor', 'k');
end

exportgraphics(fig4, ...
    fullfile(saveDir, ...
    'pairwise_nochunk_significant_lag_vs_correlation.pdf'), ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');

%% ========================================================================
%% COMMON LIMITS: CORRELATION >= 0.2
%% ========================================================================

allCorrOver = vertcat(corrOverCell{:});
allCorrOver = allCorrOver(isfinite(allCorrOver));

allLagCorrOver = vertcat(lagCorrOverCell{:});
allLagCorrOver = allLagCorrOver(isfinite(allLagCorrOver));

if isempty(allCorrOver)

    commonCorrOverXLim = [0.2 1];

else

    xMin = min(allCorrOver);
    xMax = max(allCorrOver);

    if xMin == xMax
        xPad = 0.01;
    else
        xPad = 0.05 * (xMax - xMin);
    end

    commonCorrOverXLim = ...
        [max(0.2, xMin - xPad), xMax + xPad];
end

corrOverEdges = linspace( ...
    commonCorrOverXLim(1), ...
    commonCorrOverXLim(2), ...
    51);

if isempty(allLagCorrOver)

    commonLagCorrOverXLim = [-0.2 0.2];

else

    xMin = floor(min(allLagCorrOver) * 1000) / 1000;
    xMax = ceil(max(allLagCorrOver) * 1000) / 1000;

    commonLagCorrOverXLim = ...
        [xMin - 0.001, xMax + 0.001];
end

lagCorrOverEdges = ...
    commonLagCorrOverXLim(1):0.001:commonLagCorrOverXLim(2);

if numel(lagCorrOverEdges) < 2
    lagCorrOverEdges = linspace( ...
        commonLagCorrOverXLim(1), ...
        commonLagCorrOverXLim(2), ...
        3);
end

%% ---- common histogram y-limits ----

maxCorrCount = 0;
maxLagCount = 0;

for s = 1:nSess

    vals = corrOverCell{s};

    if ~isempty(vals)

        counts = histcounts(vals, corrOverEdges);

        if ~isempty(counts)
            maxCorrCount = max(maxCorrCount, max(counts));
        end
    end

    vals = lagCorrOverCell{s};

    if ~isempty(vals)

        counts = histcounts(vals, lagCorrOverEdges);

        if ~isempty(counts)
            maxLagCount = max(maxLagCount, max(counts));
        end
    end
end

if maxCorrCount == 0
    commonCorrOverYLim = [0 1];
else
    commonCorrOverYLim = [0 1.05 * maxCorrCount];
end

if maxLagCount == 0
    commonLagCorrOverYLim = [0 1];
else
    commonLagCorrOverYLim = [0 1.05 * maxLagCount];
end

%% ========================================================================
%% FIGURE 5: PEAK CORRELATIONS >= 0.2
%% ========================================================================

fig5 = figure( ...
    'Name', 'Pairwise no-chunk peak correlations >= 0.2', ...
    'Color', 'w', ...
    'Position', [100 100 1500 650]);

tile5 = tiledlayout( ...
    1, nSess, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

title(tile5, ...
    'Pairwise no-chunk peak correlations >= 0.2', ...
    'FontSize', titleFont);

for s = 1:nSess

    ax = nexttile(tile5, s);
    hold(ax, 'on');

    vals = corrOverCell{s};

    if ~isempty(vals)

        histogram(ax, ...
            vals, ...
            'BinEdges', corrOverEdges, ...
            'FaceColor', histColor, ...
            'EdgeColor', 'none');
    end

    xlim(ax, commonCorrOverXLim);
    ylim(ax, commonCorrOverYLim);

    xlabel(ax, ...
        'Peak correlation (>= 0.2)', ...
        'FontSize', labelFont);

    ylabel(ax, ...
        'Count', ...
        'FontSize', labelFont);

    title(ax, ...
        sprintf('%s (n=%d / %d)', ...
        animalIDs(s), ...
        numel(vals), ...
        totalPairs(s)), ...
        'FontSize', titleFont);

    box(ax, 'off');
    grid(ax, 'off');

    set(ax, ...
        'FontSize', tickFont, ...
        'LineWidth', axesLineWidth, ...
        'TickDir', 'out', ...
        'TickLength', [0.025 0.025], ...
        'XColor', 'k', ...
        'YColor', 'k');
end

exportgraphics(fig5, ...
    fullfile(saveDir, ...
    'pairwise_nochunk_peak_correlations_ge_0p2.pdf'), ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');

%% ========================================================================
%% FIGURE 6: PEAK LAGS WHERE CORRELATION >= 0.2
%% ========================================================================

fig6 = figure( ...
    'Name', 'Pairwise no-chunk peak lags where correlation >= 0.2', ...
    'Color', 'w', ...
    'Position', [100 100 1500 650]);

tile6 = tiledlayout( ...
    1, nSess, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

title(tile6, ...
    'Pairwise no-chunk peak lags where peak correlation >= 0.2', ...
    'FontSize', titleFont);

for s = 1:nSess

    ax = nexttile(tile6, s);
    hold(ax, 'on');

    vals = lagCorrOverCell{s};

    if ~isempty(vals)

        histogram(ax, ...
            vals, ...
            'BinEdges', lagCorrOverEdges, ...
            'FaceColor', histColor, ...
            'EdgeColor', 'none');
    end

    xlim(ax, commonLagCorrOverXLim);
    ylim(ax, commonLagCorrOverYLim);

    xlabel(ax, ...
        'Peak lag (s)', ...
        'FontSize', labelFont);

    ylabel(ax, ...
        'Count', ...
        'FontSize', labelFont);

    title(ax, ...
        sprintf('%s (n=%d)', ...
        animalIDs(s), ...
        numel(vals)), ...
        'FontSize', titleFont);

    box(ax, 'off');
    grid(ax, 'off');

    set(ax, ...
        'FontSize', tickFont, ...
        'LineWidth', axesLineWidth, ...
        'TickDir', 'out', ...
        'TickLength', [0.025 0.025], ...
        'XColor', 'k', ...
        'YColor', 'k');
end

exportgraphics(fig6, ...
    fullfile(saveDir, ...
    'pairwise_nochunk_peak_lags_corr_ge_0p2.pdf'), ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');

%% ========================================================================
%% FIGURE 7: LAG VS CORRELATION WHERE CORRELATION >= 0.2
%% ========================================================================

fig7 = figure( ...
    'Name', 'Pairwise no-chunk lag vs correlation where correlation >= 0.2', ...
    'Color', 'w', ...
    'Position', [100 100 1500 650]);

tile7 = tiledlayout( ...
    1, nSess, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

title(tile7, ...
    'Pairwise no-chunk lag vs. correlation where peak correlation >= 0.2', ...
    'FontSize', titleFont);

for s = 1:nSess

    ax = nexttile(tile7, s);
    hold(ax, 'on');

    lagVals = lagCorrOverCell{s};
    corrVals = corrOverCell{s};

    if ~isempty(lagVals)

        scatter(ax, ...
            lagVals, ...
            corrVals, ...
            12, ...
            scatterColor, ...
            'filled');
    end

    xlim(ax, commonLagCorrOverXLim);
    ylim(ax, commonCorrOverXLim);

    xlabel(ax, ...
        'Peak lag (s)', ...
        'FontSize', labelFont);

    ylabel(ax, ...
        'Peak correlation', ...
        'FontSize', labelFont);

    title(ax, ...
        sprintf('%s (n=%d)', ...
        animalIDs(s), ...
        numel(corrVals)), ...
        'FontSize', titleFont);

    box(ax, 'off');
    grid(ax, 'off');

    set(ax, ...
        'FontSize', tickFont, ...
        'LineWidth', axesLineWidth, ...
        'TickDir', 'out', ...
        'TickLength', [0.025 0.025], ...
        'XColor', 'k', ...
        'YColor', 'k');
end

exportgraphics(fig7, ...
    fullfile(saveDir, ...
    'pairwise_nochunk_lag_vs_corr_ge_0p2.pdf'), ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');

fprintf('\nAll pairwise no-chunk paper plots saved to:\n%s\n', saveDir);

end
