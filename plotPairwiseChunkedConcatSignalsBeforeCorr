function plotPairwiseChunkedConcatSignalsBeforeCorr(dataFile, channelsToUse, intIdx, pyrIdx, doBaselineNorm, useShiftedInt)
% plots the exact concatenated interneuron and pyramidal signals before pairwise xcorr
% intIdx and pyrIdx are local indices within intCxWinCell / pyrCxWinCell
% j run for actual concat signals before corr: plotPairwiseChunkedConcatSignalsBeforeCorr("Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData\EMG_Neural_AllChannels.mat",1:4,1,1,true,false)
% j run for shift-int control signal before corr: plotPairwiseChunkedConcatSignalsBeforeCorr("Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData\EMG_Neural_AllChannels.mat",1:4,1,1,true,true)

if nargin < 2 || isempty(channelsToUse), channelsToUse = 1:4; end
if nargin < 3 || isempty(intIdx), intIdx = 1; end
if nargin < 4 || isempty(pyrIdx), pyrIdx = 1; end
if nargin < 5 || isempty(doBaselineNorm), doBaselineNorm = true; end
if nargin < 6 || isempty(useShiftedInt), useShiftedInt = false; end

S = load(dataFile, 'pyrCxWinCell', 'intCxWinCell', 'intCxWinShiftedCell', 'tAxis');

pyrEvents = [];
intEvents = [];

for ch = channelsToUse
    pyrChunk = squeeze(S.pyrCxWinCell{ch}(:, pyrIdx, :));

    if useShiftedInt
        intChunk = squeeze(S.intCxWinShiftedCell{ch}(:, intIdx, :));
    else
        intChunk = squeeze(S.intCxWinCell{ch}(:, intIdx, :));
    end

    pyrEvents = cat(1, pyrEvents, pyrChunk);
    intEvents = cat(1, intEvents, intChunk);
end

if doBaselineNorm
    pyrEvents = subtractTrialBaseline(pyrEvents, S.tAxis, -500, -450);
    intEvents = subtractTrialBaseline(intEvents, S.tAxis, -500, -450);
end

pyrConcat = pyrEvents';
intConcat = intEvents';

pyrConcat = pyrConcat(:);
intConcat = intConcat(:);

validIdx = ~isnan(pyrConcat) & ~isnan(intConcat);
pyrConcat = pyrConcat(validIdx);
intConcat = intConcat(validIdx);

maxPlotSamples = min(numel(pyrConcat), 10000);
plotIdx = 1:maxPlotSamples;

figure('Color','w','Name','Concatenated Signals Before Pairwise XCorr');
plot(plotIdx, pyrConcat(plotIdx), 'b', 'LineWidth', 1.2); hold on;
plot(plotIdx, intConcat(plotIdx), 'r', 'LineWidth', 1.2);

xlabel('Concatenated time samples');
ylabel('Firing rate');
if doBaselineNorm
    normTxt = 'baseline-subtracted';
else
    normTxt = 'raw';
end

if useShiftedInt
    shiftTxt = 'shifted interneuron vs unshifted pyramidal';
else
    shiftTxt = 'unshifted interneuron vs unshifted pyramidal';
end

title(sprintf('Before XCorr: %s, %s, Int %d vs Pyr %d', shiftTxt, normTxt, intIdx, pyrIdx));
legend({'Pyramidal','Interneuron'}, 'Location','best');
box off;

fprintf('\n===== concatenated signal check =====\n');
fprintf('intIdx = %d | pyrIdx = %d\n', intIdx, pyrIdx);
fprintf('useShiftedInt = %d | doBaselineNorm = %d\n', useShiftedInt, doBaselineNorm);
fprintf('number of valid concatenated samples = %d\n', numel(pyrConcat));
fprintf('corr at zero lag = %.6f\n', corr(intConcat, pyrConcat, 'Rows','complete'));

end

function Eout = subtractTrialBaseline(Ein, tAxis, tStart, tEnd)
[~, iStart] = min(abs(tAxis - tStart));
[~, iEnd] = min(abs(tAxis - tEnd));

if iEnd < iStart
    tmp = iStart;
    iStart = iEnd;
    iEnd = tmp;
end

baseline = mean(Ein(:, iStart:iEnd), 2, 'omitnan');
Eout = Ein - baseline;
end
