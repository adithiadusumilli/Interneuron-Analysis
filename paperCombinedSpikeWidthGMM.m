function paperCombinedSpikeWidthGMM(baseDirs, regionToPlot)
% Paper plot using the GMM saved by classifyCortexAndStriatum_getMouseDataNames.
% Histogram: D026, D020, D024, D050, D054 (D043 also included because classifications made based on this in total)
% The intersection marker is calculated from the TWO WEIGHTED curves plotted.
%
% Run:
%   paperCombinedSpikeWidthGMM([], "cortex")
%   paperCombinedSpikeWidthGMM([], "striatum")

arguments
    baseDirs = []
    regionToPlot (1,1) string = "cortex"
end
regionToPlot = lower(regionToPlot);
if ~ismember(regionToPlot, ["cortex","striatum"])
    error('regionToPlot must be "cortex" or "striatum".');
end

matFile = 'X:\David\AnalysesData\AA_classifications.mat';
saveDir = 'X:\David\AnalysesData\InterneuronAnalyses\AA Paper Plots\classifications';
if ~isfile(matFile), error('Cannot find %s',matFile); end
S = load(matFile);
wanted = ["D026","D020","D024","D050","D054"];
if ~isfield(S,'mouseIDs') || ~isfield(S,'fileWidths')
    error('AA_classifications.mat must contain mouseIDs and fileWidths.');
end
savedIDs = string(S.mouseIDs(:));
rows = zeros(1,numel(wanted));
for k = 1:numel(wanted)
    matches = find(strcmpi(savedIDs,wanted(k)));
    if numel(matches) ~= 1
        error('Expected one row for %s; found %d.',wanted(k),numel(matches));
    end
    rows(k) = matches;
end
if ~isempty(baseDirs)
    if isstring(baseDirs), baseDirs = cellstr(baseDirs); end
    if ~iscell(baseDirs), error('baseDirs must be a cell array or [].'); end
    for k = 1:numel(baseDirs)
        if contains(upper(string(baseDirs{k})),'D043')
            error('D043 is excluded from the histogram.');
        end
    end
    warning('paperCombinedSpikeWidthGMM:UnusedPaths', ...
        'baseDirs is accepted for compatibility; saved fileWidths are used.');
end

switch regionToPlot
    case "cortex"
        regionCol = 1;
        paramsField = 'gmmParamsCortex';
        lowerField = 'lowerBoundCortex';
        upperField = 'upperBoundCortex';
        wideName = 'Pyramidal fit';
        narrowName = 'Interneuron fit';
    case "striatum"
        regionCol = 2;
        paramsField = 'gmmParamsStriat';
        lowerField = 'lowerBoundStriat';
        upperField = 'upperBoundStriat';
        wideName = 'Wide-spiking fit';
        narrowName = 'Narrow-spiking fit';
end
for field = {paramsField,lowerField,upperField}
    if ~isfield(S,field{1}), error('Missing %s in %s.',field{1},matFile); end
end

widths = [];
for k = 1:numel(rows)
    w = double(S.fileWidths{rows(k),regionCol});
    widths = [widths, w(:)']; %#ok<AGROW>
end
widths = widths(isfinite(widths) & widths >= 0 & widths <= 0.8e-3);
if isempty(widths), error('No valid spike widths found.'); end
widthMs = widths * 1000;

p = S.(paramsField);
mu = double(p.mu(:)') * 1000;
sigma = double(p.sigma(:)') * 1000;
weights = double(p.weights(:)');
[mu,ord] = sort(mu);
sigma = sigma(ord);
weights = weights(ord);
if numel(mu) ~= 2 || numel(sigma) ~= 2 || numel(weights) ~= 2 || ...
        any(sigma <= 0) || any(weights <= 0)
    error('Invalid saved GMM parameters.');
end

% Weighted intersection: precisely where the two DISPLAYED component curves
% cross. The saved classification function uses an UNWEIGHTED intersection
% for reference; that is intentionally not used for this marker.
normalPdf = @(x,m,s) exp(-0.5*((x-m)./s).^2) ./ (s*sqrt(2*pi));
diffWeighted = @(x) weights(1)*normalPdf(x,mu(1),sigma(1)) - ...
                    weights(2)*normalPdf(x,mu(2),sigma(2));
scanX = linspace(mu(1),mu(2),2001);
scanY = diffWeighted(scanX);
crossings = find(scanY(1:end-1).*scanY(2:end) <= 0);
if isempty(crossings)
    error('Weighted Gaussian curves do not cross between their means.');
end
[~,best] = min(abs((scanX(crossings)+scanX(crossings+1))/2 - mean(mu)));
i = crossings(best);
intersectionMs = fzero(diffWeighted,[scanX(i),scanX(i+1)]);

lowerMs = double(S.(lowerField)) * 1000;
upperMs = double(S.(upperField)) * 1000;

fontSize = 10;
axesLineWidth = 0.5;
tickLength = [0.015 0.015];  % shorter than previous [0.025 0.025]
binWidthMs = 0.05;
xLimits = [0.15 0.8];
blue = [0 0 1];
red = [1 0 0];
lowerColor = [0 0.55 0.45]; % teal
upperColor = [0.55 0.25 0.70]; % purple

fig = figure('Color','w','Position',[100 100 740 440]);
ax = axes('Parent',fig); hold(ax,'on');
% Explicit bin edges anchor bins to zero; xlim only crops the view.
% Thin black bar edges distinguish neighboring histogram bins.
h = histogram(ax,widthMs,'BinEdges',0:binWidthMs:0.8, ...
    'FaceColor',blue,'FaceAlpha',0.55, ...
    'EdgeColor','k','LineWidth',0.5);

xMs = linspace(0,0.8,1500);
n = numel(widthMs);
yNarrow = normalPdf(xMs,mu(1),sigma(1)) * n * binWidthMs * weights(1);
yWide = normalPdf(xMs,mu(2),sigma(2)) * n * binWidthMs * weights(2);
hWide = plot(ax,xMs,yWide,'Color',blue,'LineWidth',1.5);
hNarrow = plot(ax,xMs,yNarrow,'Color',red,'LineWidth',1.5);
hIntersection = xline(ax,intersectionMs,'k--','LineWidth',1.2);
hLower = xline(ax,lowerMs,':','Color',lowerColor,'LineWidth',1.5);
hUpper = xline(ax,upperMs,'-.','Color',upperColor,'LineWidth',1.5);

xlabel(ax,'Trough-to-peak width (ms)');
ylabel(ax,'Count');
title(ax,sprintf('%s spike-width distribution (5 animals)', ...
    upperFirst(regionToPlot)),'FontWeight','normal');
legend(ax,[h,hWide,hNarrow,hIntersection,hLower,hUpper], ...
    {'Spike widths',wideName,narrowName,'Weighted curve intersection', ...
     'Lower classification bound','Upper classification bound'}, ...
    'Location','northeastoutside','Box','off','FontSize',fontSize);
xlim(ax,xLimits);
set(ax,'FontSize',fontSize,'LineWidth',axesLineWidth,'TickDir','out', ...
    'TickLength',tickLength,'Box','off','XColor','k','YColor','k','Color','w');
grid(ax,'off');
set(findall(fig,'Type','text'),'FontSize',fontSize);

if ~isfolder(saveDir), mkdir(saveDir); end
stem = sprintf('CombinedSpikeWidthGMM_%s_D026_D020_D024_D050_D054_SavedModel',char(regionToPlot));
exportgraphics(fig,fullfile(saveDir,[stem '.pdf']), ...
    'ContentType','vector','BackgroundColor','white');
savefig(fig,fullfile(saveDir,[stem '.fig']));

fprintf('\n%s widths: %d neurons from %s (D043 excluded from histogram).\n', ...
    char(regionToPlot),n,strjoin(cellstr(wanted),', '));
fprintf('Weighted plotted-curve intersection: %.4f ms\n',intersectionMs);
fprintf('Saved classification bounds: %.4f ms and %.4f ms\n',lowerMs,upperMs);
fprintf('NOTE: Saved GMM may include D043 if included during original classification.\n');
fprintf('Saved vector PDF and FIG to: %s\n',saveDir);
end

function result = upperFirst(str)
result = char(str);
result(1) = upper(result(1));
end
