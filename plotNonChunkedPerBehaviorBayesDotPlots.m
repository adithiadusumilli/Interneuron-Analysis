function plotNonChunkedPerBehaviorBayesDotPlots(bayesFile)
% One figure per animal; one subplot/tile per behavior.
% Each tile plots the two saved evidence ratios:
%   H0/H+50 and H0/H-50
% with BF=1 reference line.
%
% Title includes behavior name and nTimepoints from the original saved file.
%
% RUN:
%   plotNonChunkedPerBehaviorBayesDotPlots

arguments
    bayesFile (1,1) string = ...
        "X:\David\AnalysesData\nonchunked_xcorr_by_classifier_cortex_all6Sessions_saved_PeakLagBayesPosthoc50ms.mat"
end

if ~isfile(bayesFile)
    error('Bayes file not found:\n%s',bayesFile);
end

S = load(bayesFile,'posthocResults','summaryTable');

if ~isfield(S,'posthocResults')
    error('Bayes file does not contain posthocResults.');
end
posthocResults = S.posthocResults;

if isfield(S,'summaryTable')
    summaryTable = S.summaryTable;
elseif isfield(posthocResults,'summaryTable')
    summaryTable = posthocResults.summaryTable;
else
    error('Bayes file does not contain summaryTable.');
end

% Load original no-chunk saved file only to recover nTimepoints.
sourceResults = [];
if isfield(posthocResults,'sourceFile') && ...
        ~isempty(posthocResults.sourceFile) && ...
        isfile(posthocResults.sourceFile)

    R = load(posthocResults.sourceFile,'results');
    if isfield(R,'results')
        sourceResults = R.results;
    end
else
    warning('Could not open posthocResults.sourceFile; nTimepoints unavailable.');
end

summaryTable.Animal = string(summaryTable.Animal);
summaryTable.BehaviorName = string(summaryTable.BehaviorName);

animals = unique(summaryTable.Animal,'stable');

% Common y limits across all plots.
allBF = [ ...
    double(summaryTable.Evidence_H0_over_H50(:)); ...
    double(summaryTable.Evidence_H0_over_Hneg50(:)); ...
    1];

allBF = allBF(isfinite(allBF));

if isempty(allBF)
    error('No finite Bayes factors found.');
end

yMin = min(allBF);
yMax = max(allBF);

if yMin == yMax
    pad = max(0.1,0.1*abs(yMin));
else
    pad = 0.08*(yMax-yMin);
end

commonYLim = [max(0,yMin-pad), yMax+pad];
commonYLim(1) = min(commonYLim(1),0.9);
commonYLim(2) = max(commonYLim(2),1.1);

for a = 1:numel(animals)

    mouseID = animals(a);

    T = summaryTable(summaryTable.Animal == mouseID,:);

    [~,ord] = sort(double(T.Behavior));
    T = T(ord,:);

    nBeh = height(T);
    nCols = 4;
    nRows = ceil(nBeh/nCols);

    fig = figure( ...
        'Name',sprintf('%s Non-chunked per-behavior Bayes factors',mouseID), ...
        'Color','w', ...
        'Position',[80 80 1500 300*nRows]);

    tl = tiledlayout(fig,nRows,nCols, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    title(tl, ...
        sprintf('%s | Non-chunked population-average Bayes factors',mouseID), ...
        'FontSize',16, ...
        'Interpreter','none');

    for r = 1:nBeh

        beh = double(T.Behavior(r));
        behName = string(T.BehaviorName(r));

        bfPlus = double(T.Evidence_H0_over_H50(r));
        bfMinus = double(T.Evidence_H0_over_Hneg50(r));

        nTimepoints = getNTimepoints(sourceResults,mouseID,beh);

        nexttile;
        hold on;

        if isfinite(bfPlus) && isfinite(bfMinus)
            plot([1 2],[bfPlus bfMinus],'-', ...
                'Color',[0.75 0.75 0.75], ...
                'LineWidth',1);
        end

        scatter(1,bfPlus,90,'filled','Marker','o');
        scatter(2,bfMinus,90,'filled','Marker','^');

        yline(1,'k--','LineWidth',1.25);

        xlim([0.5 2.5]);
        ylim(commonYLim);

        xticks([1 2]);
        xticklabels({'H0/H+50','H0/H-50'});
        ylabel('Bayes factor');

        if isfinite(nTimepoints)
            title(sprintf('%s | n = %d timepoints', ...
                char(behName),round(nTimepoints)), ...
                'Interpreter','none','FontSize',11);
        else
            title(sprintf('%s | n unavailable',char(behName)), ...
                'Interpreter','none','FontSize',11);
        end

        grid on;
        box off;
        set(gca,'FontSize',10,'TickDir','out');
    end
end

end

function nTimepoints = getNTimepoints(sourceResults,mouseID,beh)

nTimepoints = NaN;

if isempty(sourceResults) || ...
        ~isfield(sourceResults,'sessions') || ...
        isempty(sourceResults.sessions)
    return;
end

sessions = sourceResults.sessions;

for s = 1:numel(sessions)

    if ~isfield(sessions(s),'mouseID') || ...
            string(sessions(s).mouseID) ~= string(mouseID)
        continue;
    end

    if ~isfield(sessions(s),'beh') || isempty(sessions(s).beh)
        return;
    end

    behVals = [sessions(s).beh.beh];
    bIdx = find(behVals == beh,1,'first');

    if isempty(bIdx)
        return;
    end

    B = sessions(s).beh(bIdx);

    if isfield(B,'nTimepoints') && ...
            ~isempty(B.nTimepoints) && ...
            isscalar(B.nTimepoints)

        nTimepoints = double(B.nTimepoints);
    end

    return;
end

end
