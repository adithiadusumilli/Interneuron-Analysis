function paperNonChunkedPerBehaviorBayesDotPlots(bayesFile)
% Non-chunked population-average Bayes factors, one figure per animal.
% D043 excluded; log y-axis; consistent paper formatting; vector PDF + FIG.
% Prints per-behavior trial counts and timepoint counts when saved.
%
% RUN: paperNonChunkedPerBehaviorBayesDotPlots

arguments
    bayesFile (1,1) string = ...
        "X:\David\AnalysesData\nonchunked_xcorr_by_classifier_cortex_all6Sessions_saved_PeakLagBayesPosthoc50ms.mat"
end

outDir = "X:\David\AnalysesData\InterneuronAnalyses\AA Paper Plots\per behavior bayes factors";
if ~isfolder(outDir), mkdir(outDir); end
if ~isfile(bayesFile), error('Bayes file not found:\n%s',bayesFile); end

S = load(bayesFile,'posthocResults','summaryTable');
if ~isfield(S,'posthocResults')
    error('Bayes file does not contain posthocResults.');
end
posthocResults = S.posthocResults;
if isfield(S,'summaryTable') && ~isempty(S.summaryTable)
    T = S.summaryTable;
elseif isfield(posthocResults,'summaryTable') && ~isempty(posthocResults.summaryTable)
    T = posthocResults.summaryTable;
else
    error('Bayes file does not contain summaryTable.');
end

required = {'Animal','Behavior','BehaviorName', ...
    'Evidence_H0_over_H50','Evidence_H0_over_Hneg50'};
for k = 1:numel(required)
    if ~ismember(required{k},T.Properties.VariableNames)
        error('summaryTable missing required variable "%s".',required{k});
    end
end
T.Animal = string(T.Animal);
T.BehaviorName = string(T.BehaviorName);
T = T(T.Animal ~= "D043",:);

% Load original results to retrieve saved nTimepoints and trial counts.
sourceResults = [];
if isfield(posthocResults,'sourceFile') && ~isempty(posthocResults.sourceFile) ...
        && isfile(string(posthocResults.sourceFile))
    R = load(string(posthocResults.sourceFile),'results');
    if isfield(R,'results'), sourceResults = R.results; end
else
    warning('BayesPlots:SourceUnavailable', ...
        'Original source file unavailable; saved trial/timepoint counts may be unavailable.');
end

n = height(T);
nTimepoints = nan(n,1);
nTrials = nan(n,1);
for i = 1:n
    [nTimepoints(i),nTrials(i)] = getBehaviorCounts(sourceResults, ...
        T.Animal(i),double(T.Behavior(i)));
end
% Prefer directly saved counts in the summary table when present.
if ismember('nTimepoints',T.Properties.VariableNames)
    v = double(T.nTimepoints);
    nTimepoints(isfinite(v)) = v(isfinite(v));
end
if ismember('nTrials',T.Properties.VariableNames)
    v = double(T.nTrials);
    nTrials(isfinite(v)) = v(isfinite(v));
end

% Print ALL behaviors, including ones without plottable evidence ratios.
outputTable = table(T.Animal,T.Behavior,T.BehaviorName,nTrials,nTimepoints, ...
    double(T.Evidence_H0_over_H50),double(T.Evidence_H0_over_Hneg50), ...
    'VariableNames',{'Animal','Behavior','BehaviorName','nTrials','nTimepoints', ...
    'Evidence_H0_over_H50','Evidence_H0_over_Hneg50'});
outputTable = sortrows(outputTable,{'Animal','Behavior'});
fprintf('\nNON-CHUNKED PER-BEHAVIOR BAYES FACTORS AND COUNTS (D043 excluded)\n');
disp(outputTable);
if any(isnan(nTrials))
    warning('BayesPlots:TrialsUnavailable', ...
        'Some nTrials values are NaN: no verified trial count was found in the saved results.');
end

plus = double(T.Evidence_H0_over_H50);
minus = double(T.Evidence_H0_over_Hneg50);
validPlus = isfinite(plus) & plus > 0;
validMinus = isfinite(minus) & minus > 0;
validRows = validPlus | validMinus;
if ~any(validRows)
    error('No positive finite Bayes factors found after excluding D043.');
end
if any((isfinite(plus) & plus <= 0) | (isfinite(minus) & minus <= 0))
    warning('BayesPlots:NonPositive', ...
        'Zero/negative Bayes factors cannot be plotted on a log axis; values remain in the printed table.');
end

allBF = [plus(validPlus);minus(validMinus);1];
lo = min(allBF); hi = max(allBF);
% Multiplicative padding, with BF=1 always visible.
logLo = log10(lo); logHi = log10(hi);
span = max(logHi-logLo,0.5);
commonYLim = 10.^[logLo-0.08*span,logHi+0.08*span];

animals = unique(T.Animal,'stable');
preferred = ["D026","D020","D024","D050","D054"];
animals = [preferred(ismember(preferred,animals)), ...
    animals(~ismember(animals,preferred)).'];
animals = animals(:);

blue = [0 0.2 0.6];
coral = [0.95 0.45 0.35];
for a = 1:numel(animals)
    mouseID = animals(a);
    idx = T.Animal == mouseID & validRows;
    A = T(idx,:);
    if isempty(A), continue; end
    [~,ord] = sort(double(A.Behavior));
    A = A(ord,:);
    originalRows = find(idx);
    originalRows = originalRows(ord);
    nBeh = height(A);
    nCols = 4;
    nRows = ceil(nBeh/nCols);
    fig = figure('Name',sprintf('%s Non-chunked per-behavior Bayes factors',mouseID), ...
        'Color','w','Position',[80 80 1500 max(300,280*nRows)]);
    tl = tiledlayout(fig,nRows,nCols,'TileSpacing','compact','Padding','compact');
    title(tl,sprintf('%s | Non-chunked population-average Bayes factors',mouseID), ...
        'FontSize',10,'FontWeight','normal','Interpreter','none');

    for r = 1:nBeh
        i = originalRows(r);
        bfPlus = plus(i); bfMinus = minus(i);
        ax = nexttile(tl);
        hold(ax,'on');
        if validPlus(i) && validMinus(i)
            plot(ax,[1 2],[bfPlus bfMinus],'-', ...
                'Color',[0.78 0.78 0.78],'LineWidth',0.8,'HandleVisibility','off');
        end
        % NaN makes an invalid point invisible while keeping legend handles.
        if ~validPlus(i), bfPlus = NaN; end
        if ~validMinus(i), bfMinus = NaN; end
        hPlus = scatter(ax,1,bfPlus,45,blue,'filled','Marker','o');
        hMinus = scatter(ax,2,bfMinus,45,coral,'filled','Marker','o');
        hRef = yline(ax,1,'k--','LineWidth',1);
        set(ax,'YScale','log','YLim',commonYLim, ...
            'XLim',[0.5 2.5],'XTick',[1 2], ...
            'XTickLabel',{'H0/H+50','H0/H-50'}, ...
            'FontSize',10,'TickDir','out','TickLength',[0.015 0.015], ...
            'LineWidth',0.5,'XColor','k','YColor','k','Color','w', ...
            'Box','off');
        grid(ax,'off');
        ylabel(ax,'Bayes factor (log scale)');
        if isfinite(nTimepoints(i))
            ttl = sprintf('%s | n = %d timepoints', ...
                char(A.BehaviorName(r)),round(nTimepoints(i)));
        else
            ttl = sprintf('%s | n unavailable',char(A.BehaviorName(r)));
        end
        title(ax,ttl,'Interpreter','none','FontSize',10,'FontWeight','normal');
        legend(ax,[hPlus,hMinus,hRef], ...
            {'H0 / H+50','H0 / H-50','BF = 1 (equal evidence)'}, ...
            'Location','best','Box','off','FontSize',8);
    end

    pdfFile = fullfile(outDir,mouseID + "_NonChunkedPerBehaviorBayesFactors.pdf");
    figFile = fullfile(outDir,mouseID + "_NonChunkedPerBehaviorBayesFactors.fig");
    exportgraphics(fig,pdfFile,'ContentType','vector','BackgroundColor','white');
    savefig(fig,figFile);
    fprintf('Saved %s\nSaved %s\n',pdfFile,figFile);
end
end

function [nTimepoints,nTrials] = getBehaviorCounts(sourceResults,mouseID,beh)
nTimepoints = NaN;
nTrials = NaN;
if isempty(sourceResults) || ~isstruct(sourceResults) || ...
        ~isfield(sourceResults,'sessions') || isempty(sourceResults.sessions)
    return;
end
sessions = sourceResults.sessions;
for s = 1:numel(sessions)
    if iscell(sessions), R = sessions{s}; else, R = sessions(s); end
    if isempty(R) || ~isstruct(R), continue; end
    if isfield(R,'mouseID'), id = string(R.mouseID);
    elseif isfield(R,'animalID'), id = string(R.animalID);
    else, continue; end
    if id ~= string(mouseID), continue; end
    if ~isfield(R,'beh') || isempty(R.beh), return; end
    Bs = R.beh;
    if iscell(Bs), Bs = [Bs{:}]; end
    if ~isstruct(Bs) || ~isfield(Bs,'beh'), return; end
    bIdx = find([Bs.beh] == beh,1,'first');
    if isempty(bIdx), return; end
    B = Bs(bIdx);
    if isfield(B,'nTimepoints') && isnumeric(B.nTimepoints) && ...
            isscalar(B.nTimepoints) && isfinite(B.nTimepoints)
        nTimepoints = double(B.nTimepoints);
    end
    % Only use fields explicitly identifying a trial count or trial IDs.
    if isfield(B,'nTrials') && isnumeric(B.nTrials) && ...
            isscalar(B.nTrials) && isfinite(B.nTrials)
        nTrials = double(B.nTrials);
    elseif isfield(B,'numTrials') && isnumeric(B.numTrials) && ...
            isscalar(B.numTrials) && isfinite(B.numTrials)
        nTrials = double(B.numTrials);
    elseif isfield(B,'trialIDs') && ~isempty(B.trialIDs)
        nTrials = numel(unique(B.trialIDs));
    elseif isfield(B,'trialIndices') && ~isempty(B.trialIndices)
        nTrials = numel(unique(B.trialIndices));
    end
    return;
end
end
