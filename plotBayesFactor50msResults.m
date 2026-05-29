function plotBayesFactor50msResults(combinedFile)
% plots BF-style 50 ms shift results from combined trial-avg chunked file
% j run: plotBayesFactor50msResults("C:\Users\mirilab\Documents\GlobusTransfer\concatCrossCorr_trialavg_chunked_popavg_ALL_SESSIONS.mat")

if nargin < 1 || isempty(combinedFile)
    combinedFile = "C:\Users\mirilab\Documents\GlobusTransfer\concatCrossCorr_trialavg_chunked_popavg_ALL_SESSIONS.mat";
end

S = load(combinedFile, 'allSessions');
sessions = S.allSessions.sessions;
nSess = numel(sessions);

animalIDs = strings(nSess,1);
bfVals = nan(nSess,1);
pOrig = nan(nSess,1);
pShift = nan(nSess,1);
lagOrig = nan(nSess,1);
lagShift = nan(nSess,1);

for i = 1:nSess
    animalIDs(i) = string(sessions(i).sessionTag);

    if isfield(sessions(i),'bf_H0_over_H1')
        bfVals(i) = sessions(i).bf_H0_over_H1;
        pOrig(i) = sessions(i).pValOriginal;
        pShift(i) = sessions(i).pValShift50;
        lagOrig(i) = sessions(i).peakLagSecOriginal;
        lagShift(i) = sessions(i).peakLagSecShift50;
    end
end

%% print table
T = table(animalIDs, bfVals, pOrig, pShift, lagOrig, lagShift, ...
    'VariableNames', {'Animal','BF_H0_over_H1','pValOriginal','pValShift50','PeakLagOriginal','PeakLagShift50'});

disp(T);

%% plot bayes factor
figure('Color','w','Name','Bayes Factor 50 ms Shift Results');
bar(bfVals);
hold on;
yline(1,'k--','LineWidth',1.5);

xticks(1:nSess);
xticklabels(animalIDs);
ylabel('BF-style ratio H0/H1');
xlabel('Animal');
title('Evidence Ratio: Original Timing vs 50 ms Interneuron Shift');
box off;
set(gca,'FontSize',14,'LineWidth',1,'TickDir','out');

%% plot p-values
figure('Color','w','Name','Bayes Factor Input p-values');
bar([pOrig pShift]);
xticks(1:nSess);
xticklabels(animalIDs);
ylabel('Permutation p-value');
xlabel('Animal');
legend({'Original timing','50 ms shifted interneuron'}, 'Location','best');
title('Permutation p-values Used in BF-style Ratio');
box off;
set(gca,'FontSize',14,'LineWidth',1,'TickDir','out');

%% plot peak lags
figure('Color','w','Name','Original vs 50 ms Shifted Peak Lags');
bar([lagOrig lagShift]);
hold on;
yline(0,'k:','LineWidth',1.5);

xticks(1:nSess);
xticklabels(animalIDs);
ylabel('Peak lag (s)');
xlabel('Animal');
legend({'Original timing','50 ms shifted interneuron'}, 'Location','best');
title('Peak Lags Before and After 50 ms Interneuron Shift');
box off;
set(gca,'FontSize',14,'LineWidth',1,'TickDir','out');

end
