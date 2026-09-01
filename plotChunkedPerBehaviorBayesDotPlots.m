function plotChunkedPerBehaviorBayesDotPlots(bayesFile)
% Plot ONLY the saved Bayes-factor dot plots for the CHUNKED
% trial-averaged population-average per-behavior analysis.
%
% Makes ONE FIGURE PER ANIMAL.
% Within each animal figure, there is ONE TILE PER BEHAVIOR.
%
% Each behavior tile contains two dots:
%   x = H0/H+50 and H0/H-50
%   y = saved Bayes factor / evidence ratio
%
% A horizontal dashed line at BF = 1 marks equal evidence.
%
% Tile title includes:
%   behavior name
%   n = number of behavior events/trials used in the chunked analysis
%
% This function DOES NOT recompute Bayes factors.
% It loads the already-saved chunked per-behavior Bayes output and follows
% posthocResults.sourceFile back to the combined chunked behavior file only
% to retrieve the saved trial/event counts.
%
% DEFAULT BAYES FILE:
%   C:\Users\mirilab\Documents\GlobusTransfer\
%   concatCrossCorrPerCanonicalBehavior_classifier_chunked_trialavg_popavg_ALL_SESSIONS_PeakLagBayesPosthoc50ms.mat
%
% RUN:
%   plotChunkedPerBehaviorBayesDotPlots

arguments
    bayesFile (1,1) string = ...
        "C:\Users\mirilab\Documents\GlobusTransfer\concatCrossCorrPerCanonicalBehavior_classifier_chunked_trialavg_popavg_ALL_SESSIONS_PeakLagBayesPosthoc50ms.mat"
end

%% ---------------- load saved Bayes output ----------------

if ~isfile(bayesFile)
    error('Chunked Bayes file not found:\n%s',bayesFile);
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
    summaryTable = [];
end

%% ---------------- load original chunked combined file for event counts ----------------

sourceAllSessions = [];

if isfield(posthocResults,'sourceFile') && ...
        ~isempty(posthocResults.sourceFile) && ...
        isfile(posthocResults.sourceFile)

    R = load(posthocResults.sourceFile,'allSessions');

    if isfield(R,'allSessions')
        sourceAllSessions = R.allSessions;
    end
else
    warning(['Could not open posthocResults.sourceFile. Event/trial counts ' ...
             'will be shown as unavailable.']);
end

%% ---------------- canonical behavior names ----------------

canonicalNames = { ...
    'climbdown', ...
    'climbup', ...
    'eating', ...
    'grooming', ...
    'jumpdown', ...
    'jumping', ...
    'rearing', ...
    'still', ...
    'walkflat', ...
    'walkgrid'};

%% ---------------- extract Bayes results into arrays ----------------

[animals,behaviors,bfH50,bfHneg50] = ...
    extractChunkedBayesValues(posthocResults,summaryTable);

if isempty(animals)
    error('Could not extract animal/behavior Bayes-factor values.');
end

uniqueAnimals = unique(animals,'stable');
uniqueBehaviors = unique(behaviors(isfinite(behaviors)),'stable');

% Prefer canonical 1:10 ordering.
orderedBehaviors = [];
for b = 1:10
    if any(uniqueBehaviors == b)
        orderedBehaviors(end+1) = b; %#ok<AGROW>
    end
end
remaining = uniqueBehaviors(~ismember(uniqueBehaviors,orderedBehaviors));
orderedBehaviors = [orderedBehaviors(:); remaining(:)]';

%% ---------------- common y limits across all animals/behaviors ----------------

allBF = [bfH50(:); bfHneg50(:); 1];
allBF = allBF(isfinite(allBF));

if isempty(allBF)
    error('No finite Bayes-factor values found.');
end

yMin = min(allBF);
yMax = max(allBF);

if yMin == yMax
    pad = max(0.1,0.10*abs(yMin));
else
    pad = 0.08*(yMax-yMin);
end

commonYLim = [max(0,yMin-pad), yMax+pad];
commonYLim(1) = min(commonYLim(1),0.9);
commonYLim(2) = max(commonYLim(2),1.1);

%% ========================================================================
%% one figure per animal
%% ========================================================================

for a = 1:numel(uniqueAnimals)

    mouseID = uniqueAnimals(a);

    nBeh = numel(orderedBehaviors);
    nCols = 5;
    nRows = ceil(nBeh/nCols);

    fig = figure( ...
        'Name',sprintf('%s Chunked per-behavior Bayes factors',mouseID), ...
        'Color','w', ...
        'Position',[60 80 1700 330*nRows]);

    tl = tiledlayout(fig,nRows,nCols, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    title(tl, ...
        sprintf('%s | Chunked trial-averaged population Bayes factors',mouseID), ...
        'FontSize',16, ...
        'Interpreter','none');

    for bIdx = 1:nBeh

        beh = orderedBehaviors(bIdx);

        rowIdx = find(animals == mouseID & behaviors == beh,1,'first');

        if isempty(rowIdx)
            thisBFPlus = NaN;
            thisBFMinus = NaN;
        else
            thisBFPlus = bfH50(rowIdx);
            thisBFMinus = bfHneg50(rowIdx);
        end

        if beh >= 1 && beh <= 10
            behName = canonicalNames{beh};
        else
            behName = sprintf('behavior %d',beh);
        end

        nEvents = getChunkedEventCount(sourceAllSessions,mouseID,beh);

        nexttile;
        hold on;

        if isfinite(thisBFPlus) && isfinite(thisBFMinus)
            plot([1 2],[thisBFPlus thisBFMinus],'-', ...
                'Color',[0.75 0.75 0.75], ...
                'LineWidth',1);
        end

        scatter(1,thisBFPlus,90,'filled','Marker','o');
        scatter(2,thisBFMinus,90,'filled','Marker','^');

        yline(1,'k--','LineWidth',1.25);

        xlim([0.5 2.5]);
        ylim(commonYLim);

        xticks([1 2]);
        xticklabels({'H0/H+50','H0/H-50'});
        ylabel('Bayes factor');

        if isfinite(nEvents)
            title(sprintf('%s | n = %d events', ...
                behName,round(nEvents)), ...
                'Interpreter','none', ...
                'FontSize',11);
        else
            title(sprintf('%s | n unavailable',behName), ...
                'Interpreter','none', ...
                'FontSize',11);
        end

        grid on;
        box off;
        set(gca,'FontSize',10,'TickDir','out');
    end
end

end

%% ========================================================================
%% helper: extract saved chunked Bayes values
%% ========================================================================

function [animals,behaviors,bfH50,bfHneg50] = ...
    extractChunkedBayesValues(posthocResults,summaryTable)

animals = strings(0,1);
behaviors = zeros(0,1);
bfH50 = zeros(0,1);
bfHneg50 = zeros(0,1);

% First preference: summary table, if present.
if ~isempty(summaryTable) && istable(summaryTable)

    vars = summaryTable.Properties.VariableNames;

    animalVar = '';
    behaviorVar = '';
    plusVar = '';
    minusVar = '';

    if ismember('Animal',vars), animalVar = 'Animal'; end
    if ismember('Behavior',vars), behaviorVar = 'Behavior'; end
    if ismember('Evidence_H0_over_H50',vars), plusVar = 'Evidence_H0_over_H50'; end
    if ismember('Evidence_H0_over_Hneg50',vars), minusVar = 'Evidence_H0_over_Hneg50'; end

    if ~isempty(animalVar) && ~isempty(behaviorVar) && ...
            ~isempty(plusVar) && ~isempty(minusVar)

        animals = string(summaryTable.(animalVar));
        behaviors = double(summaryTable.(behaviorVar));
        bfH50 = double(summaryTable.(plusVar));
        bfHneg50 = double(summaryTable.(minusVar));
        return;
    end
end

% Fallback: posthocResults.behaviorResults{bIdx}(animalIdx)
if isfield(posthocResults,'behaviorResults') && ...
        ~isempty(posthocResults.behaviorResults)

    BR = posthocResults.behaviorResults;

    for bIdx = 1:numel(BR)

        X = BR{bIdx};

        if isempty(X)
            continue;
        end

        for a = 1:numel(X)

            R = X(a);

            if isfield(R,'mouseID') && ~isempty(R.mouseID)
                mouseID = string(R.mouseID);
            elseif isfield(R,'animalID') && ~isempty(R.animalID)
                mouseID = string(R.animalID);
            else
                mouseID = "Animal_" + a;
            end

            if isfield(R,'beh') && ~isempty(R.beh)
                beh = double(R.beh);
            elseif isfield(posthocResults,'behaviors') && ...
                    numel(posthocResults.behaviors) >= bIdx
                beh = double(posthocResults.behaviors(bIdx));
            else
                beh = bIdx;
            end

            plus = NaN;
            minus = NaN;

            if isfield(R,'evidenceRatio_H0_over_H50')
                plus = double(R.evidenceRatio_H0_over_H50);
            elseif isfield(R,'Evidence_H0_over_H50')
                plus = double(R.Evidence_H0_over_H50);
            end

            if isfield(R,'evidenceRatio_H0_over_Hneg50')
                minus = double(R.evidenceRatio_H0_over_Hneg50);
            elseif isfield(R,'Evidence_H0_over_Hneg50')
                minus = double(R.Evidence_H0_over_Hneg50);
            end

            animals(end+1,1) = mouseID; %#ok<AGROW>
            behaviors(end+1,1) = beh; %#ok<AGROW>
            bfH50(end+1,1) = plus; %#ok<AGROW>
            bfHneg50(end+1,1) = minus; %#ok<AGROW>
        end
    end
end

end

%% ========================================================================
%% helper: retrieve real behavior event/trial count
%% ========================================================================

function nEvents = getChunkedEventCount(allSessions,mouseID,beh)

nEvents = NaN;

if isempty(allSessions) || ...
        ~isfield(allSessions,'behaviorResults') || ...
        isempty(allSessions.behaviorResults)

    return;
end

% Find behavior cell by saved .beh value.
for bIdx = 1:numel(allSessions.behaviorResults)

    B = allSessions.behaviorResults{bIdx};

    if isempty(B) || ~isfield(B,'beh') || isempty(B.beh)
        continue;
    end

    savedBeh = B.beh;
    if numel(savedBeh) > 1
        savedBeh = savedBeh(1);
    end

    if double(savedBeh) ~= double(beh)
        continue;
    end

    % Find animal row.
    animalIdx = [];

    if isfield(allSessions,'mouseIDs') && ~isempty(allSessions.mouseIDs)

        ids = string(allSessions.mouseIDs);
        animalIdx = find(ids == string(mouseID),1,'first');
    end

    if isempty(animalIdx)
        return;
    end

    % Preferred count: real_nTrialsUsed if present.
    if isfield(B,'real_nTrialsUsed') && ...
            numel(B.real_nTrialsUsed) >= animalIdx && ...
            isfinite(B.real_nTrialsUsed(animalIdx))

        nEvents = double(B.real_nTrialsUsed(animalIdx));
        return;
    end

    % Fallback: total number of real trials/events for that behavior.
    if isfield(B,'real_nTrialsThisBehavior') && ...
            numel(B.real_nTrialsThisBehavior) >= animalIdx && ...
            isfinite(B.real_nTrialsThisBehavior(animalIdx))

        nEvents = double(B.real_nTrialsThisBehavior(animalIdx));
        return;
    end

    return;
end

end
