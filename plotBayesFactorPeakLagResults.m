function plotBayesFactorPeakLagResults()

combinedFile = "C:\Users\mirilab\Documents\GlobusTransfer\concatCrossCorr_trialavg_chunked_popavg_ALL_SESSIONS.mat";

S = load(combinedFile,'allSessions');

sessions = S.allSessions.sessions;
nSess = numel(sessions);

animalIDs = strings(nSess,1);
BF = nan(nSess,1);
pH0 = nan(nSess,1);
pH50 = nan(nSess,1);
realLag = nan(nSess,1);

for i = 1:nSess

    animalIDs(i) = string(sessions(i).sessionTag);

    BF(i)      = sessions(i).evidenceRatio_H0_over_H50;
    pH0(i)     = sessions(i).pValH0;
    pH50(i)    = sessions(i).pValH50;
    realLag(i) = sessions(i).bayes_realPeakLagSec;

end

%% summary table

T = table( ...
    animalIDs, ...
    realLag, ...
    pH0, ...
    pH50, ...
    BF, ...
    'VariableNames', ...
    {'Animal','RealPeakLag_s','pValH0','pValH50','BF_H0_over_H50'});

disp(T)

%% ---------------- BF plot ----------------

figure('Color','w');
bar(BF)

hold on
yline(1,'k--','LineWidth',1.5)

xticks(1:nSess)
xticklabels(animalIDs)

ylabel('Evidence Ratio H0 / H50')
xlabel('Animal')

title('Bayes-style Evidence Ratio')

box off
set(gca,'FontSize',14)

%% ---------------- p-value comparison ----------------

figure('Color','w');

bar([pH0 pH50])

xticks(1:nSess)
xticklabels(animalIDs)

ylabel('Permutation p-value')
xlabel('Animal')

legend({'H0 (0 ms)','H50 (+50 ms lead)'}, ...
    'Location','best')

title('Observed Peak Lag Under Competing Models')

box off
set(gca,'FontSize',14)

%% ---------------- real lag plot ----------------

figure('Color','w');

bar(realLag)

hold on
yline(0,'k--')
yline(0.050,'r--')

xticks(1:nSess)
xticklabels(animalIDs)

ylabel('Real Peak Lag (s)')
xlabel('Animal')

title('Observed Peak Lag Relative to H0 and H50')

legend({'Real Peak Lag','0 ms','50 ms'}, ...
    'Location','best')

box off
set(gca,'FontSize',14)

end
