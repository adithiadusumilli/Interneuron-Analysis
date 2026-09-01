function plotPairwiseLagImbalanceResults_FromGlobusTransfer( ...
    noChunkFile, chunkFile)
% Plot saved pairwise lag-imbalance results OUTSIDE QUEST.
%
% Uses the FINAL DAVID-UPDATE outputs directly.
%
% DEFAULT INPUTS:
%
% C:\Users\mirilab\Documents\GlobusTransfer\
%   pairwiseNoChunkLagImbalanceBayes50ms_StoreyCorrThresh_DAVID_UPDATE.mat
%
% C:\Users\mirilab\Documents\GlobusTransfer\
%   pairwiseChunkedLagImbalanceBayes50ms_StoreyCorrThresh_DAVID_UPDATE.mat
%
% RUN:
%   plotPairwiseLagImbalanceResults_FromGlobusTransfer
%
% Figures made separately for NO-CHUNK and CHUNKED:
%
%   1. Actual lag imbalance with 95% H0 permutation bounds
%   2. Per-animal H0 permutation distributions
%   3. Combined H0 permutation distribution
%   4. H0 / H+50 / H-50 distributions per animal
%
% Lag imbalance:
%
%   (# positive lags - # negative lags)
%   ---------------------------------
%   (# positive lags + # negative lags)
%
% Exact zero lags are excluded.

arguments
    noChunkFile (1,1) string = ...
        "C:\Users\mirilab\Documents\GlobusTransfer\pairwiseNoChunkLagImbalanceBayes50ms_StoreyCorrThresh_DAVID_UPDATE.mat"

    chunkFile (1,1) string = ...
        "C:\Users\mirilab\Documents\GlobusTransfer\pairwiseChunkedLagImbalanceBayes50ms_StoreyCorrThresh_DAVID_UPDATE.mat"
end


%% ========================================================================
%% SETTINGS
%% ========================================================================

preferredAnimalOrder = [ ...
    "D026", ...
    "D020", ...
    "D024", ...
    "D043", ...
    "D050", ...
    "D054"];

actualColor = [0.95 0.45 0.35];
permHistColor = [0.3 0.6 0.8];

h0Color   = [0.30 0.60 0.80];
h50Color  = [0.85 0.45 0.25];
hm50Color = [0.45 0.70 0.40];

ciColor = [0.20 0.20 0.20];

legendFont = 16;
titleFont  = 18;
labelFont  = 18;
tickFont   = 16;


%% ========================================================================
%% LOAD
%% ========================================================================

noChunkResults = loadResults(noChunkFile, "no-chunk");
chunkResults   = loadResults(chunkFile,   "chunked");


%% ========================================================================
%% PLOT BOTH ANALYSES
%% ========================================================================

plotOneAnalysis( ...
    noChunkResults, ...
    preferredAnimalOrder, ...
    "NO-CHUNK", ...
    actualColor, ...
    permHistColor, ...
    h0Color, ...
    h50Color, ...
    hm50Color, ...
    ciColor, ...
    legendFont, ...
    titleFont, ...
    labelFont, ...
    tickFont);

plotOneAnalysis( ...
    chunkResults, ...
    preferredAnimalOrder, ...
    "CHUNKED TRIAL-AVERAGED", ...
    actualColor, ...
    permHistColor, ...
    h0Color, ...
    h50Color, ...
    hm50Color, ...
    ciColor, ...
    legendFont, ...
    titleFont, ...
    labelFont, ...
    tickFont);

end


%% ========================================================================
%% LOAD RESULTS
%% ========================================================================

function results = loadResults(fileName, analysisLabel)

if ~isfile(fileName)
    error('%s file not found:\n%s', analysisLabel, fileName);
end

S = load(fileName, 'results');

if ~isfield(S, 'results')
    error('%s file does not contain variable "results".', analysisLabel);
end

results = S.results;

if ~isfield(results, 'sessions') || isempty(results.sessions)
    error('%s results.sessions is missing or empty.', analysisLabel);
end

end


%% ========================================================================
%% PLOT ONE ANALYSIS
%% ========================================================================

function plotOneAnalysis( ...
    results, preferredAnimalOrder, analysisLabel, ...
    actualColor, permHistColor, ...
    h0Color, h50Color, hm50Color, ciColor, ...
    legendFont, titleFont, labelFont, tickFont)


%% ------------------------------------------------------------------------
%% EXTRACT SAVED RESULTS
%% ------------------------------------------------------------------------

sessions = results.sessions;
nSessRaw = numel(sessions);

animalIDsRaw = strings(nSessRaw,1);

actualRaw = nan(nSessRaw,1);
nSigRaw   = nan(nSessRaw,1);

H0Raw   = cell(nSessRaw,1);
H50Raw  = cell(nSessRaw,1);
Hm50Raw = cell(nSessRaw,1);

ciH0Raw   = nan(nSessRaw,2);
ciH50Raw  = nan(nSessRaw,2);
ciHm50Raw = nan(nSessRaw,2);

bfPlusRaw  = nan(nSessRaw,1);
bfMinusRaw = nan(nSessRaw,1);


for s = 1:nSessRaw

    if iscell(sessions)
        R = sessions{s};
    else
        R = sessions(s);
    end

    if isempty(R)
        animalIDsRaw(s) = "Animal_" + s;
        continue;
    end


    %% animal ID

    if isfield(R,'animalID') && ~isempty(R.animalID)

        animalIDsRaw(s) = string(R.animalID);

    elseif isfield(R,'mouseID') && ~isempty(R.mouseID)

        animalIDsRaw(s) = string(R.mouseID);

    elseif isfield(R,'baseDir') && ~isempty(R.baseDir)

        hit = regexp(char(R.baseDir),'D\d+','match','once');

        if isempty(hit)
            animalIDsRaw(s) = "Animal_" + s;
        else
            animalIDsRaw(s) = string(hit);
        end

    else

        animalIDsRaw(s) = "Animal_" + s;

    end


    %% actual lag imbalance

    if isfield(R,'actualLagImbalance')
        actualRaw(s) = double(R.actualLagImbalance);
    end


    %% number of actual significant int-pyr pairs

    if isfield(R,'nActualSigIntPyr')
        nSigRaw(s) = double(R.nActualSigIntPyr);
    end


    %% permutation distributions
    %
    % Prefer the already-filtered valid distributions.

    if isfield(R,'validNullLagImbalanceH0')
        H0Raw{s} = double(R.validNullLagImbalanceH0(:));
    elseif isfield(R,'nullLagImbalanceH0')
        x = double(R.nullLagImbalanceH0(:));
        H0Raw{s} = x(isfinite(x));
    else
        H0Raw{s} = [];
    end


    if isfield(R,'validNullLagImbalanceH50')
        H50Raw{s} = double(R.validNullLagImbalanceH50(:));
    elseif isfield(R,'nullLagImbalanceH50')
        x = double(R.nullLagImbalanceH50(:));
        H50Raw{s} = x(isfinite(x));
    else
        H50Raw{s} = [];
    end


    if isfield(R,'validNullLagImbalanceHneg50')
        Hm50Raw{s} = double(R.validNullLagImbalanceHneg50(:));
    elseif isfield(R,'nullLagImbalanceHneg50')
        x = double(R.nullLagImbalanceHneg50(:));
        Hm50Raw{s} = x(isfinite(x));
    else
        Hm50Raw{s} = [];
    end


    %% saved 95% permutation intervals

    if isfield(R,'nullCIH0') && numel(R.nullCIH0) == 2
        ciH0Raw(s,:) = double(R.nullCIH0(:)');
    elseif ~isempty(H0Raw{s})
        ciH0Raw(s,:) = prctile(H0Raw{s}, [2.5 97.5]);
    end


    if isfield(R,'nullCIH50') && numel(R.nullCIH50) == 2
        ciH50Raw(s,:) = double(R.nullCIH50(:)');
    elseif ~isempty(H50Raw{s})
        ciH50Raw(s,:) = prctile(H50Raw{s}, [2.5 97.5]);
    end


    if isfield(R,'nullCIHneg50') && numel(R.nullCIHneg50) == 2
        ciHm50Raw(s,:) = double(R.nullCIHneg50(:)');
    elseif ~isempty(Hm50Raw{s})
        ciHm50Raw(s,:) = prctile(Hm50Raw{s}, [2.5 97.5]);
    end


    %% evidence ratios

    if isfield(R,'evidenceRatio_H0_over_H50')
        bfPlusRaw(s) = double(R.evidenceRatio_H0_over_H50);
    end

    if isfield(R,'evidenceRatio_H0_over_Hneg50')
        bfMinusRaw(s) = double(R.evidenceRatio_H0_over_Hneg50);
    end

end


%% ------------------------------------------------------------------------
%% PUT ANIMALS IN CANONICAL ORDER
%% ------------------------------------------------------------------------

animalOrder = strings(0,1);

for i = 1:numel(preferredAnimalOrder)

    if any(animalIDsRaw == preferredAnimalOrder(i))
        animalOrder(end+1,1) = preferredAnimalOrder(i); %#ok<AGROW>
    end

end

remaining = animalIDsRaw(~ismember(animalIDsRaw, animalOrder));
animalOrder = [animalOrder; remaining];

nSess = numel(animalOrder);


actual = nan(nSess,1);
nSig   = nan(nSess,1);

H0   = cell(nSess,1);
H50  = cell(nSess,1);
Hm50 = cell(nSess,1);

ciH0   = nan(nSess,2);
ciH50  = nan(nSess,2);
ciHm50 = nan(nSess,2);

bfPlus  = nan(nSess,1);
bfMinus = nan(nSess,1);


for s = 1:nSess

    idx = find(animalIDsRaw == animalOrder(s), 1, 'first');

    if isempty(idx)
        continue;
    end

    actual(s) = actualRaw(idx);
    nSig(s)   = nSigRaw(idx);

    H0{s}   = H0Raw{idx};
    H50{s}  = H50Raw{idx};
    Hm50{s} = Hm50Raw{idx};

    ciH0(s,:)   = ciH0Raw(idx,:);
    ciH50(s,:)  = ciH50Raw(idx,:);
    ciHm50(s,:) = ciHm50Raw(idx,:);

    bfPlus(s)  = bfPlusRaw(idx);
    bfMinus(s) = bfMinusRaw(idx);

end


%% ========================================================================
%% FIGURE 1
%% ACTUAL LAG IMBALANCE + H0 PERMUTATION BOUNDS
%% ========================================================================

figure( ...
    'Name', sprintf('%s Pairwise Lag Imbalance Summary', analysisLabel), ...
    'Color','w');

hold on;

x = (1:nSess)';


% 95% permutation interval
for s = 1:nSess

    if all(isfinite(ciH0(s,:)))

        line( ...
            [x(s) x(s)], ...
            ciH0(s,:), ...
            'Color',[0.60 0.60 0.60], ...
            'LineWidth',3);

    end

end


% actual values
scatter( ...
    x, ...
    actual, ...
    90, ...
    actualColor, ...
    'filled');


% zero = balanced positive/negative lag counts
yline(0,'k:','LineWidth',1.5);


xlim([0.5 nSess+0.5]);

xticks(1:nSess);
xticklabels(cellstr(animalOrder));

xlabel('Animal','FontSize',labelFont);
ylabel('Lag Imbalance','FontSize',labelFont);

title( ...
    sprintf('%s Pairwise Lag Imbalance with 95%% H0 Permutation Bounds', ...
    analysisLabel), ...
    'FontSize',titleFont);

box off;
grid on;

set(gca, ...
    'FontSize',tickFont, ...
    'LineWidth',1, ...
    'TickDir','out');


%% add exact actual values above dots

yl = ylim;
yr = diff(yl);

for s = 1:nSess

    if isfinite(actual(s))

        text( ...
            x(s), ...
            actual(s) + 0.025*yr, ...
            sprintf('%.3f',actual(s)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'FontSize',12);

    end

end


%% ========================================================================
%% DETERMINE COMMON HISTOGRAM LIMITS
%% ========================================================================

allH0 = [];

for s = 1:nSess

    vals = H0{s};

    if ~isempty(vals)
        allH0 = [allH0; vals(:)]; %#ok<AGROW>
    end

end

allH0 = allH0(isfinite(allH0));

actualFinite = actual(isfinite(actual));

allVals = [allH0(:); actualFinite(:)];

if isempty(allVals)

    commonXLim = [-1 1];

else

    xMin = min(allVals);
    xMax = max(allVals);

    if xMin == xMax
        pad = 0.05;
    else
        pad = 0.05*(xMax-xMin);
    end

    commonXLim = [xMin-pad xMax+pad];

end

% Lag imbalance is mathematically bounded by -1 and +1.
commonXLim(1) = max(commonXLim(1), -1);
commonXLim(2) = min(commonXLim(2),  1);

nBins = 24;

commonEdges = linspace( ...
    commonXLim(1), ...
    commonXLim(2), ...
    nBins+1);


%% ========================================================================
%% FIGURE 2
%% PER-ANIMAL H0 PERMUTATION DISTRIBUTIONS
%% ========================================================================

figure( ...
    'Name', sprintf('%s H0 Lag-Imbalance Permutations',analysisLabel), ...
    'Color','w');

tl = tiledlayout( ...
    1,nSess, ...
    'TileSpacing','compact', ...
    'Padding','compact');

title( ...
    tl, ...
    sprintf('%s: Actual Lag Imbalance vs. H0 Permutation Distribution', ...
    analysisLabel), ...
    'FontSize',titleFont);


legendSet = false;

for s = 1:nSess

    ax = nexttile;
    hold(ax,'on');

    vals = H0{s};

    vals = vals(isfinite(vals));

    if isempty(vals)

        title( ...
            sprintf('%s (no valid H0 draws)',animalOrder(s)), ...
            'FontSize',titleFont);

        axis off;
        continue;

    end


    hHist = histogram( ...
        vals, ...
        'BinEdges',commonEdges, ...
        'FaceColor',permHistColor, ...
        'EdgeColor','none');


    bounds = ciH0(s,:);

    hLow = gobjects(1);
    hHigh = gobjects(1);

    if all(isfinite(bounds))

        hLow = xline( ...
            bounds(1), ...
            '--', ...
            'Color',ciColor, ...
            'LineWidth',1.5);

        hHigh = xline( ...
            bounds(2), ...
            '--', ...
            'Color',ciColor, ...
            'LineWidth',1.5);

    end


    hActual = gobjects(1);

    if isfinite(actual(s))

        hActual = xline( ...
            actual(s), ...
            '-', ...
            'Color',actualColor, ...
            'LineWidth',2);

    end


    xlim(commonXLim);

    xlabel('Lag Imbalance','FontSize',labelFont);
    ylabel('Count','FontSize',labelFont);

    if isfinite(nSig(s))

        title( ...
            sprintf('%s | n_{sig} = %d', ...
            animalOrder(s),round(nSig(s))), ...
            'FontSize',titleFont);

    else

        title(animalOrder(s),'FontSize',titleFont);

    end

    box off;

    set(gca, ...
        'FontSize',tickFont, ...
        'LineWidth',1, ...
        'TickDir','out');


    if ~legendSet && isgraphics(hActual)

        handles = hHist;
        labels = {'H0 Permuted Lag Imbalance'};

        handles(end+1) = hActual;
        labels{end+1} = 'Actual Lag Imbalance';

        if isgraphics(hLow)
            handles(end+1) = hLow;
            labels{end+1} = '2.5% H0 Bound';
        end

        if isgraphics(hHigh)
            handles(end+1) = hHigh;
            labels{end+1} = '97.5% H0 Bound';
        end

        lgd = legend( ...
            handles, ...
            labels, ...
            'Orientation','horizontal');

        lgd.Layout.Tile = 'south';
        lgd.FontSize = legendFont;
        lgd.Box = 'off';

        legendSet = true;

    end

end


%% ========================================================================
%% FIGURE 3
%% COMBINED H0 PERMUTATION DISTRIBUTION
%% ========================================================================

if ~isempty(allH0)

    figure( ...
        'Name',sprintf('%s Combined H0 Permutation Distribution',analysisLabel), ...
        'Color','w');

    hold on;


    hHist = histogram( ...
        allH0, ...
        'BinEdges',commonEdges, ...
        'FaceColor',permHistColor, ...
        'EdgeColor','none');


    pooledCI = prctile(allH0,[2.5 97.5]);

    hLow = xline( ...
        pooledCI(1), ...
        '--', ...
        'Color',ciColor, ...
        'LineWidth',1.5);

    hHigh = xline( ...
        pooledCI(2), ...
        '--', ...
        'Color',ciColor, ...
        'LineWidth',1.5);


    co = lines(nSess);

    actualLines = gobjects(0,1);
    actualLabels = {};

    for s = 1:nSess

        if isfinite(actual(s))

            actualLines(end+1,1) = xline( ...
                actual(s), ...
                '-', ...
                'Color',co(s,:), ...
                'LineWidth',1.8); %#ok<AGROW>

            actualLabels{end+1} = sprintf( ...
                '%s Actual = %.3f', ...
                animalOrder(s), ...
                actual(s)); %#ok<AGROW>

        end

    end


    xlim(commonXLim);

    xlabel('Lag Imbalance','FontSize',labelFont);
    ylabel('Count','FontSize',labelFont);

    title( ...
        sprintf('%s: All Animals Combined H0 Permutation Distribution', ...
        analysisLabel), ...
        'FontSize',titleFont);

    box off;

    set(gca, ...
        'FontSize',tickFont, ...
        'LineWidth',1, ...
        'TickDir','out');


    legHandles = [hHist hLow hHigh actualLines(:)'];

    legLabels = { ...
        'H0 Permuted Lag Imbalance', ...
        sprintf('2.5%% pooled = %.3f',pooledCI(1)), ...
        sprintf('97.5%% pooled = %.3f',pooledCI(2))};

    legLabels = [legLabels actualLabels];

    legend( ...
        legHandles, ...
        legLabels, ...
        'Location','best', ...
        'Box','off', ...
        'FontSize',legendFont);

end


%% ========================================================================
%% FIGURE 4
%% H0 / H+50 / H-50 DISTRIBUTIONS PER ANIMAL
%% ========================================================================

figure( ...
    'Name',sprintf('%s Three-Model Permutation Distributions',analysisLabel), ...
    'Color','w');

tl3 = tiledlayout( ...
    1,nSess, ...
    'TileSpacing','compact', ...
    'Padding','compact');

title( ...
    tl3, ...
    sprintf('%s: H0, H+50, and H-50 Lag-Imbalance Distributions', ...
    analysisLabel), ...
    'FontSize',titleFont);


% Determine common limits using all three models.
allThree = [];

for s = 1:nSess

    allThree = [ ...
        allThree; ...
        H0{s}(:); ...
        H50{s}(:); ...
        Hm50{s}(:)]; %#ok<AGROW>

end

allThree = allThree(isfinite(allThree));

if ~isempty(allThree)

    xMin = min([allThree(:); actualFinite(:)]);
    xMax = max([allThree(:); actualFinite(:)]);

    if xMin == xMax
        pad = 0.05;
    else
        pad = 0.05*(xMax-xMin);
    end

    modelXLim = [max(-1,xMin-pad) min(1,xMax+pad)];

else

    modelXLim = [-1 1];

end

modelEdges = linspace(modelXLim(1),modelXLim(2),25);


legendMade = false;

for s = 1:nSess

    ax = nexttile;
    hold(ax,'on');


    v0 = H0{s};
    vp = H50{s};
    vm = Hm50{s};

    v0 = v0(isfinite(v0));
    vp = vp(isfinite(vp));
    vm = vm(isfinite(vm));


    if isempty(v0) && isempty(vp) && isempty(vm)

        title( ...
            sprintf('%s (no valid draws)',animalOrder(s)), ...
            'FontSize',titleFont);

        axis off;
        continue;

    end


    h0 = histogram( ...
        v0, ...
        'BinEdges',modelEdges, ...
        'Normalization','probability', ...
        'DisplayStyle','stairs', ...
        'EdgeColor',h0Color, ...
        'LineWidth',2);


    hp = histogram( ...
        vp, ...
        'BinEdges',modelEdges, ...
        'Normalization','probability', ...
        'DisplayStyle','stairs', ...
        'EdgeColor',h50Color, ...
        'LineWidth',2);


    hm = histogram( ...
        vm, ...
        'BinEdges',modelEdges, ...
        'Normalization','probability', ...
        'DisplayStyle','stairs', ...
        'EdgeColor',hm50Color, ...
        'LineWidth',2);


    hActual = gobjects(1);

    if isfinite(actual(s))

        hActual = xline( ...
            actual(s), ...
            '-', ...
            'Color',actualColor, ...
            'LineWidth',2);

    end


    xlim(modelXLim);

    xlabel('Lag Imbalance','FontSize',labelFont);
    ylabel('Probability','FontSize',labelFont);


    if isfinite(nSig(s))

        title( ...
            sprintf('%s | n_{sig} = %d', ...
            animalOrder(s),round(nSig(s))), ...
            'FontSize',titleFont);

    else

        title(animalOrder(s),'FontSize',titleFont);

    end


    box off;

    set(gca, ...
        'FontSize',tickFont, ...
        'LineWidth',1, ...
        'TickDir','out');


    if ~legendMade && isgraphics(hActual)

        lgd = legend( ...
            [h0 hp hm hActual], ...
            {'H0','H+50','H-50','Actual Lag Imbalance'}, ...
            'Orientation','horizontal');

        lgd.Layout.Tile = 'south';
        lgd.FontSize = legendFont;
        lgd.Box = 'off';

        legendMade = true;

    end

end


%% ========================================================================
%% PRINT TABLE
%% ========================================================================

fprintf('\n\n');
fprintf('============================================================\n');
fprintf('%s PAIRWISE LAG-IMBALANCE RESULTS\n',analysisLabel);
fprintf('============================================================\n');


T = table( ...
    animalOrder, ...
    actual, ...
    ciH0(:,1), ...
    ciH0(:,2), ...
    nSig, ...
    bfPlus, ...
    bfMinus, ...
    'VariableNames',{ ...
        'Animal', ...
        'ActualLagImbalance', ...
        'H0_2p5', ...
        'H0_97p5', ...
        'NumSignificantIntPyrPairs', ...
        'Evidence_H0_over_Hplus50', ...
        'Evidence_H0_over_Hminus50'});

disp(T);

end
