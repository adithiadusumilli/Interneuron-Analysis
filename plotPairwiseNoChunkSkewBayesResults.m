function plotPairwiseNoChunkSkewBayesResults()

bayesFile = "C:\Users\mirilab\Documents\GlobusTransfer\pairwiseNoChunkSkewBayes50ms_StoreyCorrThresh.mat";

S = load(bayesFile,"results");

sessions = S.results.sessions;

sessions = sessions(~cellfun(@isempty,sessions));

nSess = numel(sessions);

animalIDs = strings(nSess,1);

BF    = nan(nSess,1);
pH0   = nan(nSess,1);
pH50  = nan(nSess,1);

realSkew = nan(nSess,1);

nSigPairs = nan(nSess,1);

for i = 1:nSess

    R = sessions{i};

    animalIDs(i) = string(R.animalID);

    BF(i) = R.evidenceRatio_H0_over_H50;

    pH0(i) = R.pValH0;
    pH50(i) = R.pValH50;

    realSkew(i) = R.actual.skew;

    nSigPairs(i) = R.actual.nSigFDR;

end

%% ---------------- summary table ----------------

T = table( ...
    animalIDs,...
    nSigPairs,...
    realSkew,...
    pH0,...
    pH50,...
    BF,...
    'VariableNames',...
    {'Animal',...
     'SignificantPairs',...
     'RealSkew',...
     'pValH0',...
     'pValH50',...
     'BF_H0_over_H50'});

disp(T)

%% ===============================================================
%% Bayes factor
%% ===============================================================

figure('Color','w');

bar(BF)

hold on

yline(1,'k--','LineWidth',1.5)

xticks(1:nSess)
xticklabels(animalIDs)

ylabel('Evidence Ratio H0 / H50')

xlabel('Animal')

title('Pairwise No-Chunk Skew Bayes Factor')

box off

set(gca,'FontSize',14)

%% ===============================================================
%% permutation p-values
%% ===============================================================

figure('Color','w');

bar([pH0 pH50])

xticks(1:nSess)
xticklabels(animalIDs)

ylabel('Tail Probability')

xlabel('Animal')

legend({'H0 (0 ms)','H50 (50 ms lead)'},...
    'Location','best')

title('Observed Skew Under Competing Models')

box off

set(gca,'FontSize',14)

%% ===============================================================
%% real skew
%% ===============================================================

figure('Color','w');

bar(realSkew)

hold on

yline(0,'k--')

xticks(1:nSess)
xticklabels(animalIDs)

ylabel('Skew')

xlabel('Animal')

title('Observed Skew of Significant Pairwise Lags')

box off

set(gca,'FontSize',14)

%% ===============================================================
%% number of significant pairs
%% ===============================================================

figure('Color','w');

bar(nSigPairs)

xticks(1:nSess)
xticklabels(animalIDs)

ylabel('# Significant Int-Pyr Pairs')

xlabel('Animal')

title('Storey + Correlation Threshold')

box off

set(gca,'FontSize',14)

end
