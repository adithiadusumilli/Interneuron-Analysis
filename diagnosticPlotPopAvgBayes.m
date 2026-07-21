function diagnosticPlotPopAvgBayes(posthocResults)

sessions = posthocResults.sessions;
nSess = numel(sessions);

%% ============================================================
%% Figure 1: overlay histograms
%% ============================================================

figure('Color','w','Name','Population Peak Lag Models');
tiledlayout(2,ceil(nSess/2),'TileSpacing','compact','Padding','compact');

for i = 1:nSess

    S = sessions(i);

    nexttile
    hold on

    allVals = [S.peakLagsH0(:);
               S.peakLagsH50(:);
               S.peakLagsHneg50(:)];

    edges = linspace(min(allVals),max(allVals),20);

    histogram(S.peakLagsH0,...
        edges,...
        'Normalization','probability',...
        'FaceColor',[0 0.45 0.74],...
        'FaceAlpha',0.35,...
        'EdgeColor','none');

    histogram(S.peakLagsH50,...
        edges,...
        'Normalization','probability',...
        'FaceColor',[0.85 0.33 0.10],...
        'FaceAlpha',0.35,...
        'EdgeColor','none');

    histogram(S.peakLagsHneg50,...
        edges,...
        'Normalization','probability',...
        'FaceColor',[0.47 0.67 0.19],...
        'FaceAlpha',0.35,...
        'EdgeColor','none');

    xline(S.realPeakLagSec,'k','LineWidth',2)

    title(S.mouseID)
    xlabel('Peak lag (s)')
    ylabel('Probability')

    legend({'H0','H+50','H-50','Real'},...
        'Location','best',...
        'Box','off')

    box off

end

%% ============================================================
%% Figure 2: boxplots
%% ============================================================

figure('Color','w','Name','Distribution summaries');

tiledlayout(2,ceil(nSess/2),'TileSpacing','compact','Padding','compact');

for i = 1:nSess

    S = sessions(i);

    nexttile
    hold on

    vals = [S.peakLagsH0(:);
            S.peakLagsH50(:);
            S.peakLagsHneg50(:)];

    groups = [ones(numel(S.peakLagsH0),1);
              2*ones(numel(S.peakLagsH50),1);
              3*ones(numel(S.peakLagsHneg50),1)];

    boxplot(vals,groups,...
        'Labels',{'H0','H+50','H-50'})

    scatter(1,S.realPeakLagSec,80,'k','filled')

    title(S.mouseID)

    ylabel('Peak lag (s)')

end

%% ============================================================
%% Table
%% ============================================================

T = table;

for i = 1:nSess

    S = sessions(i);

    T.Animal(i) = string(S.mouseID);

    T.RealLag(i) = S.realPeakLagSec;

    T.H0Median(i) = median(S.peakLagsH0);

    T.H50Median(i) = median(S.peakLagsH50);

    T.Hneg50Median(i) = median(S.peakLagsHneg50);

    T.DistanceToH0(i) = abs(S.realPeakLagSec-T.H0Median(i));

    T.DistanceToH50(i) = abs(S.realPeakLagSec-T.H50Median(i));

    T.DistanceToHneg50(i) = abs(S.realPeakLagSec-T.Hneg50Median(i));

    T.H0overH50(i) = S.evidenceRatio_H0_over_H50;

    T.H0overHneg50(i) = S.evidenceRatio_H0_over_Hneg50;

end

disp(' ')
disp('================ Diagnostic Table ================')
disp(T)

%% ============================================================
%% Means/stds
%% ============================================================

fprintf('\nDistribution summaries\n\n')

for i = 1:nSess

    S = sessions(i);

    fprintf('%s\n',S.mouseID)

    fprintf(' H0     mean %.4f   std %.4f\n',...
        mean(S.peakLagsH0),std(S.peakLagsH0))

    fprintf(' H+50   mean %.4f   std %.4f\n',...
        mean(S.peakLagsH50),std(S.peakLagsH50))

    fprintf(' H-50   mean %.4f   std %.4f\n\n',...
        mean(S.peakLagsHneg50),std(S.peakLagsHneg50))

end

end
