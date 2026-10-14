function paperCombinedSpikeWidthGMM(baseDirs, regionToPlot)
% Paper plot of the EXISTING global GMM from
% classifyCortexAndStriatum_getMouseDataNames (saved results in AA_Classifications.mat)

% Run:
%   plotCombinedSpikeWidthGMM([], "cortex")
%   plotCombinedSpikeWidthGMM([], "striatum")

arguments
    baseDirs = []
    regionToPlot (1,1) string = "cortex"
end

regionToPlot = lower(regionToPlot);
if ~ismember(regionToPlot, ["cortex", "striatum"])
    error('regionToPlot must be "cortex" or "striatum".');
end

matFile = 'X:\David\AnalysesData\AA_classifications.mat';
if ~isfile(matFile)
    error('Classification file not found: %s', matFile);
end
S = load(matFile);
required = {'mouseIDs','fileWidths'};
for k = 1:numel(required)
    if ~isfield(S, required{k})
        error('Missing %s in %s. Rerun classifyCortexAndStriatum_getMouseDataNames.', required{k}, matFile);
    end
end

wanted = ["D026", "D020", "D024", "D050", "D054"];
savedIDs = string(S.mouseIDs(:));
selectedRows = zeros(1,numel(wanted));
for k = 1:numel(wanted)
    hits = find(strcmpi(savedIDs, wanted(k)));
    if numel(hits) ~= 1
        error('Expected exactly one entry for %s in saved mouseIDs; found %d.', wanted(k), numel(hits));
    end
    selectedRows(k) = hits;
end

if ~isempty(baseDirs)
    if isstring(baseDirs), baseDirs = cellstr(baseDirs); end
    if ~iscell(baseDirs)
        error('baseDirs must be a cell array, string array, or [].');
    end
    for k = 1:numel(baseDirs)
        thisPath = string(baseDirs{k});
        if contains(upper(thisPath), 'D043')
            error('D043 is excluded. Remove it from baseDirs or pass [].');
        end
    end
    warning('baseDirs is retained for compatibility; widths are read from the saved classification file.');
end

switch regionToPlot
    case "cortex"
        regionCol = 1;
        paramsField = 'gmmParamsCortex';
        intersectionField = 'intersectionPointCortex';
        lowerField = 'lowerBoundCortex';
        upperField = 'upperBoundCortex';
        wideName = 'Pyramidal fit';
        narrowName = 'Interneuron fit';
    case "striatum"
        regionCol = 2;
        paramsField = 'gmmParamsStriat';
        intersectionField = 'intersectionPointStriat';
        lowerField = 'lowerBoundStriat';
        upperField = 'upperBoundStriat';
        wideName = 'Wide-spiking fit';
        narrowName = 'Narrow-spiking fit';
end
for f = {paramsField, intersectionField, lowerField, upperField}
    if ~isfield(S, f{1})
        error('Missing %s in classification file.', f{1});
    end
end

% Exactly the width values used by the classification function, in seconds.
widths = [];
for k = 1:numel(selectedRows)
    w = double(S.fileWidths{selectedRows(k),regionCol});
    widths = [widths, w(:)']; %#ok<AGROW>
end
% Classification function fitted its GMM using widths <= 0.8 ms.
widths = widths(isfinite(widths) & widths >= 0 & widths <= 0.8e-3);
if isempty(widths)
    error('No valid widths found for the five selected animals.');
end

p = S.(paramsField);
mu = double(p.mu(:)');
sigma = double(p.sigma(:)');
weights = double(p.weights(:)');
[mu, idx] = sort(mu);
sigma = sigma(idx);
weights = weights(idx);
if numel(mu) ~= 2 || any(sigma <= 0)
    error('Saved GMM parameters are invalid.');
end

intersectionMs = 1e3 * S.(intersectionField);
lowerMs = 1e3 * S.(lowerField);
upperMs = 1e3 * S.(upperField);
widthMs = widths * 1e3;
binWidthMs = 0.05; % same as classification function

fig = figure('Color','w','Position',[100 100 740 440]);
ax = axes(fig); hold(ax,'on');
h = histogram(ax, widthMs, 'BinWidth',binWidthMs, ...
    'FaceColor',[0 0 1], 'FaceAlpha',0.55, 'EdgeColor','none');

xMs = linspace(0,0.8,1000);
% Saved model curves are weighted by the mixture proportions, as in the
% classification function's pooled-distribution plot.
normalPdf = @(x,m,s) exp(-0.5*((x-m)./s).^2)./(s*sqrt(2*pi));
n = numel(widthMs);
yNarrow = normalPdf(xMs, mu(1)*1e3, sigma(1)*1e3) * n * h.BinWidth * weights(1);
yWide = normalPdf(xMs, mu(2)*1e3, sigma(2)*1e3) * n * h.BinWidth * weights(2);
hWide = plot(ax,xMs,yWide,'Color',[0 0 1],'LineWidth',1.5);
hNarrow = plot(ax,xMs,yNarrow,'Color',[1 0 0],'LineWidth',1.5);
hIntersection = xline(ax,intersectionMs,'k--','LineWidth',1);
hLower = xline(ax,lowerMs,'Color',[0.25 0.25 0.25],'LineStyle',':','LineWidth',1.2);
hUpper = xline(ax,upperMs,'Color',[0.25 0.25 0.25],'LineStyle','-.','LineWidth',1.2);

xlabel(ax,'Trough-to-peak width (ms)');
ylabel(ax,'Count');
title(ax,sprintf('%s spike-width distribution (5 animals)',upperFirst(regionToPlot)), ...
    'FontWeight','normal');
legend(ax,[h,hWide,hNarrow,hIntersection,hLower,hUpper], ...
    {'Spike widths',char(wideName),char(narrowName), ...
     'Saved GMM intersection','Lower classification bound','Upper classification bound'}, ...
    'Location','northeast','Box','off');
xlim(ax,[0 0.8]);
set(ax,'FontSize',10,'LineWidth',0.5,'TickDir','out', ...
    'TickLength',[0.025 0.025],'Box','off','XColor','k','YColor','k','Color','w');
grid(ax,'off');
set(findall(fig,'Type','text'),'FontSize',10);

saveDir = 'X:\David\AnalysesData\InterneuronAnalyses\AA Paper Plots\spike width GMM';
if ~exist(saveDir,'dir'), mkdir(saveDir); end
stem = sprintf('CombinedSpikeWidthGMM_%s_D026_D020_D024_D050_D054_SavedModel',char(regionToPlot));
exportgraphics(fig,fullfile(saveDir,[stem '.pdf']), ...
    'ContentType','vector','BackgroundColor','white');
savefig(fig,fullfile(saveDir,[stem '.fig']));

fprintf('\n%s: displayed widths from %s (D043 excluded)\n',regionToPlot,strjoin(wanted,', '));
fprintf('Number of plotted widths (<=0.8 ms): %d\n',n);
fprintf('Saved model intersection: %.4f ms\n',intersectionMs);
fprintf('Saved classification boundaries: %.4f to %.4f ms\n',lowerMs,upperMs);
fprintf('NOTE: saved GMM parameters may include D043 if it was in the original classification run.\n');
fprintf('Saved PDF and FIG to: %s\n',saveDir);
end

function result = upperFirst(str)
result = char(str);
result(1) = upper(result(1));
end
