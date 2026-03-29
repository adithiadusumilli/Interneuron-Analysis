function plotSignificantPairwiseNoChunk_QCorrThresh(combinedMatFile, alpha, nNullDraws, corrThresh)

arguments
    combinedMatFile (1,1) string
    alpha (1,1) double = 0.05
    nNullDraws (1,1) double = 100
    corrThresh (1,1) double = 0.02
end

if exist('mafdr', 'file') ~= 2
    error('mafdr not available');
end

S = load(combinedMatFile, 'allSessions');
sessions = S.allSessions.sessions;
numSessions = numel(sessions);

summaryActualSkew = nan(numSessions,1);
summaryNullCI = nan(numSessions,2);
summaryAnimalID = cell(1,numSessions);

rng(0);

for sess = 1:numSessions

    peakCorrs = sessions(sess).peakCorrMatAll;
    peakLags  = sessions(sess).peakLagMatAll;
    nullXC    = sessions(sess).nullCorrMatAllShifts;

    nInt = sessions(sess).nInt;
    nPyr = sessions(sess).nPyr;
    nAll = sessions(sess).nAll;

    animalID = regexp(sessions(sess).baseDir, 'D\d+', 'match', 'once');
    if isempty(animalID)
        animalID = sprintf('Session%d', sess);
    end
    summaryAnimalID{sess} = animalID;

    [rows, cols] = getIntPyrPairs(nInt, nPyr);
    actual = computePairStatsQCorr(peakCorrs, peakLags, nullXC, rows, cols, alpha, corrThresh);

    %% ================= HEATMAP =================
    sigMat = nan(nInt, nPyr);
    for k = 1:numel(actual.rows)
        r = actual.rows(k);
        c = actual.cols(k) - nInt;
        if actual.sigMask(k)
            sigMat(r,c) = actual.lagVals(k);
        end
    end

    figure('Color','w');
    hImg = imagesc(sigMat);
    set(hImg,'AlphaData',~isnan(sigMat));
    set(gca,'Color','k');
    axis xy;
    colormap(parula);
    cb = colorbar;
    ylabel(cb,'Peak Lag (Seconds)');
    xlabel('Pyramidal Neuron Index');
    ylabel('Interneuron Index');

    title(sprintf('%s: Significant Pairwise Peak Lag Map (%d / %d Pairs)', ...
        animalID, actual.nSig, actual.nPairsNominal));

    box off; set(gca,'TickDir','out');

    %% ================= LAG HIST =================
    figure('Color','w');
    if ~isempty(actual.sigLagVec)
        edges = makeLagEdges(actual.sigLagVec);
        histogram(actual.sigLagVec,'BinEdges',edges,'EdgeColor','none');
    else
        histogram(nan);
    end
    hold on
    xline(0,'k--','LineWidth',1.5);

    xlabel('Peak Lag (Seconds)');
    ylabel('Pair Count');

    title(sprintf('%s: Distribution Of Peak Lags (%d Significant Pairs)', ...
        animalID, actual.nSig));

    grid on

    %% ================= SCATTER =================
    if ~isempty(actual.sigLagVec)
        figure('Color','w');
        scatterhist(actual.sigLagVec(:), actual.sigCorrVec(:), ...
            'Direction','out','Marker','.');

        hold on
        yline(corrThresh,'r--','LineWidth',1.5);

        xlabel('Peak Lag (Seconds)');
        ylabel('Peak Correlation');

        title(sprintf('%s: Peak Lag Versus Peak Correlation (%d Significant Pairs)', ...
            animalID, actual.nSig));

        legend({'Pairs','Correlation Threshold'},'Location','best');
    end

    %% ================= SKEW =================
    [poolRows, poolCols] = getAllUpperPairs(nAll);
    nullSkews = nan(nNullDraws,1);

    for r = 1:nNullDraws
        idx = randperm(numel(poolRows), actual.nPairsNominal);
        nullDraw = computePairStatsQCorr(peakCorrs, peakLags, nullXC, ...
            poolRows(idx), poolCols(idx), alpha, corrThresh);
        nullSkews(r) = nullDraw.skew;
    end

    validNull = nullSkews(~isnan(nullSkews));

    if isempty(validNull)
        nullCI = [NaN NaN];
    else
        nullCI = prctile(validNull,[2.5 97.5]);
    end

    summaryActualSkew(sess) = actual.skew;
    summaryNullCI(sess,:) = nullCI;

    figure('Color','w');
    histogram(validNull,30,'EdgeColor','none'); hold on
    xline(actual.skew,'r','LineWidth',2);
    xline(nullCI(1),'k--','LineWidth',1.5);
    xline(nullCI(2),'k--','LineWidth',1.5);

    xlabel('Skew (Mean Minus Median Over Standard Deviation)');
    ylabel('Count');

    title(sprintf('%s: Skew Relative To Null Distribution', animalID));

    legend({'Null Distribution','Actual Skew','Null 95% Bounds'}, ...
        'Location','best');

    grid on

end

%% ================= SUMMARY =================
figure('Color','w');
tiledlayout(1,numSessions,'TileSpacing','compact','Padding','compact');

for s = 1:numSessions
    nexttile; hold on

    xNull = 0.95;
    xActual = 1.05;

    if ~any(isnan(summaryNullCI(s,:)))
        line([xNull xNull], summaryNullCI(s,:), ...
            'Color',[0.5 0.5 0.5],'LineWidth',3);
    end

    plot(xActual, summaryActualSkew(s),'ko','MarkerFaceColor','k');

    yline(0,'k:');

    xlim([0.7 1.3]);
    xticks(1);
    xticklabels({summaryAnimalID{s}});

    ylabel('Skew');

    title(summaryAnimalID{s});

    box off; grid on
end

title(t,'Summary Of Skew Relative To Null Distribution Across Sessions');

end
