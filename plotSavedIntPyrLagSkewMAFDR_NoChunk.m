function plotSavedIntPyrLagSkewMAFDR_NoChunk(resultsFile)

arguments
    resultsFile (1,1) string
end

S = load(resultsFile, 'results');
R = S.results;

nSess = numel(R.sessions);

summaryActualSkew = nan(nSess,1);
summaryNullCI = nan(nSess,2);
summaryAnimalID = cell(1,nSess);
summaryNSig = nan(nSess,1);
summaryNTotal = nan(nSess,1);

for s = 1:nSess
    sess = R.sessions{s};

    animalID = sess.animalID;
    summaryAnimalID{s} = animalID;

    actual = sess.actual;

    nInt = sess.nInt;
    nPyr = sess.nPyr;

    % ================= HEATMAP =================
    sigMat = nan(nInt, nPyr);

    for k = 1:numel(actual.rows)
        r = actual.rows(k);
        c = actual.cols(k) - nInt;

        if c >= 1 && c <= nPyr && actual.sigFDRMask(k)
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
    ylabel(cb,'Peak Lag (Seconds)','FontSize',14);

    xlabel('Pyramidal Neuron Index','FontSize',14);
    ylabel('Interneuron Index','FontSize',14);

    title(sprintf('%s Significant Pairwise Peak Lag Map (%d / %d Pairs)', ...
        animalID, actual.nSigFDR, actual.nPairsNominal), ...
        'FontSize',16,'FontWeight','bold');

    set(gca,'FontSize',13,'TickDir','out');
    box off;

    % ================= LAG HIST =================
    figure('Color','w');
    if ~isempty(actual.sigLagVec)
        edges = makeLagEdges(actual.sigLagVec);
        histogram(actual.sigLagVec,'BinEdges',edges,'EdgeColor','none');
    else
        histogram(nan);
    end
    hold on

    hZero = xline(0,'k--','LineWidth',1.5);

    xlabel('Peak Lag (Seconds)','FontSize',14);
    ylabel('Pair Count','FontSize',14);

    title(sprintf('%s Distribution Of Peak Lags (%d / %d Significant Pairs)', ...
        animalID, actual.nSigFDR, actual.nPairsNominal), ...
        'FontSize',16,'FontWeight','bold');

    legend(hZero,'Zero Lag Reference','Location','best','FontSize',12);

    set(gca,'FontSize',13);
    grid on

    % ================= SCATTER =================
    if ~isempty(actual.sigLagVec)
        sigCorrVec = actual.realVals(actual.sigFDRMask);
        sigCorrVec = sigCorrVec(~isnan(sigCorrVec) & isfinite(sigCorrVec));

        if numel(sigCorrVec) == numel(actual.sigLagVec)
            figure('Color','w');
            scatterhist(actual.sigLagVec(:), sigCorrVec(:), ...
                'Direction','out','Marker','.');

            hold on
            hThresh = yline(R.corrThresh,'r--','LineWidth',1.5);

            xlabel('Peak Lag (Seconds)','FontSize',14);
            ylabel('Peak Correlation','FontSize',14);

            title(sprintf('%s Peak Lag Versus Peak Correlation (%d / %d Significant Pairs)', ...
                animalID, actual.nSigFDR, actual.nPairsNominal), ...
                'FontSize',16,'FontWeight','bold');

            legend(hThresh,'Correlation Threshold','FontSize',12,'Location','best');

            set(gca,'FontSize',13);
        end
    end

    % ================= SKEW =================
    validNullSkews = sess.validNullSkews;

    if ~isempty(validNullSkews)
        nullCI = prctile(validNullSkews,[2.5 97.5]);
    else
        nullCI = [NaN NaN];
    end

    summaryActualSkew(s) = actual.skew;
    summaryNullCI(s,:) = nullCI;
    summaryNSig(s) = actual.nSigFDR;
    summaryNTotal(s) = actual.nPairsNominal;

    figure('Color','w');

    if ~isempty(validNullSkews)
        histogram(validNullSkews,30,'EdgeColor','none'); hold on
        hActual = xline(actual.skew,'r','LineWidth',2);
        hCI = xline(nullCI(1),'k--','LineWidth',1.5);
        xline(nullCI(2),'k--','LineWidth',1.5);

        legend([hActual hCI], ...
            {'Actual Skew','95% Null Skew Interval'}, ...
            'FontSize',12,'Location','best');
    else
        histogram(nan);
    end

    xlabel('Skew (Mean Minus Median Over Standard Deviation)','FontSize',14);
    ylabel('Count','FontSize',14);

    title(sprintf('%s Skew Relative To Null Distribution',animalID), ...
        'FontSize',16,'FontWeight','bold');

    set(gca,'FontSize',13);
    grid on
end

% ================= SUMMARY =================
figure('Color','w');
t = tiledlayout(1,nSess,'TileSpacing','compact','Padding','compact');

for s = 1:nSess
    nexttile; hold on

    xNull = 0.95;
    xActual = 1.05;

    hCI = gobjects(1);
    if ~any(isnan(summaryNullCI(s,:)))
        hCI = line([xNull xNull], summaryNullCI(s,:), ...
            'Color',[0.5 0.5 0.5],'LineWidth',4);
    end

    hDot = plot(xActual, summaryActualSkew(s), ...
        'ko','MarkerFaceColor','k','MarkerSize',9);

    yline(0,'k:','LineWidth',1.5);

    xlim([0.7 1.3]);
    xticks(1);
    xticklabels({summaryAnimalID{s}});

    ylabel('Skew','FontSize',16);

    title(sprintf('%s (%d / %d)', ...
        summaryAnimalID{s}, summaryNSig(s), summaryNTotal(s)), ...
        'FontSize',16,'FontWeight','bold');

    set(gca,'FontSize',15);
    box off
    grid on

    if s == 1 && isgraphics(hCI)
        lgd = legend([hDot hCI], ...
            {'Actual Skew','95% Null Skew Interval'}, ...
            'Location','southoutside','FontSize',16);
        lgd.Layout.Tile = 'south';
    end
end

title(t,'Summary Of Skew Relative To Null Distribution Across Sessions', ...
    'FontSize',18,'FontWeight','bold');

end
