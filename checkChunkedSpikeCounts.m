function spikeStats = checkChunkedSpikeCounts(dataFile, channelsToUse)
% checks whether shifted EMG-aligned neural windows have different spike/activity counts than unshifted windows
% re-extracts shifted counts for all shifts from validTransitionsNeurShiftedCell
% assuming tAxis is in ms
% j run: checkChunkedSpikeCounts("Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData\EMG_Neural_AllChannels.mat",1:4)

if nargin < 2 || isempty(channelsToUse)
    channelsToUse = 1:4;
end

baseDir = fileparts(dataFile);

S = load(dataFile, ...
    'pyrCxWinCell', ...
    'intCxWinCell', ...
    'validTransitionsNeurShiftedCell', ...
    'tAxis');

dtSec = median(diff(S.tAxis)) / 1000;
tOffsets = round(S.tAxis(:)');

% load full cortex firing rates
F = load(fullfile(baseDir, 'NeuralFiringRates1msBins10msGauss.mat'), 'cortexFRs', 'cortexInds');

cortexFRs = F.cortexFRs;
cortexInds = F.cortexInds(:);

% load classifications
conslidatedDataFoler = 'X:\David\AnalysesData';
C = load(fullfile(conslidatedDataFoler, 'AA_classifications.mat'), 'classifications');
classifications = C.classifications;

animalFolders = {
    'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData', ...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData', ...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData', ...
    'X:\David\ArenaRecordings\D043-020425-ArenaRecording\ProcessedData'
};

matchRow = find(contains(string(animalFolders), string(baseDir)), 1);

if isempty(matchRow)
    error('could not match baseDir to animalFolders for AA_classifications.');
end

cortexLabelsAll = classifications{matchRow, 1};
regionClass = cortexLabelsAll(cortexInds);

pyrFRs = cortexFRs(regionClass == 0, :);
intFRs = cortexFRs(regionClass == 1, :);

% real counts from saved unshifted windows
pyrReal = [];
intReal = [];

for ch = channelsToUse
    pyrReal = cat(1, pyrReal, sum(S.pyrCxWinCell{ch}, 3, 'omitnan') * dtSec);
    intReal = cat(1, intReal, sum(S.intCxWinCell{ch}, 3, 'omitnan') * dtSec);
end

% validTransitionsNeurShiftedCell stores:
% column 1 = unshifted
% columns 2:end = shifted controls
shiftCols = 2:size(S.validTransitionsNeurShiftedCell, 2);
nShifts = numel(shiftCols);

pyrShiftMeanByShift = nan(nShifts, 1);
intShiftMeanByShift = nan(nShifts, 1);

pyrShiftCountsByShift = cell(nShifts, 1);
intShiftCountsByShift = cell(nShifts, 1);

for iShift = 1:nShifts
    shiftCol = shiftCols(iShift);

    pyrShiftThis = [];
    intShiftThis = [];

    for ch = channelsToUse
        centers = S.validTransitionsNeurShiftedCell{ch, shiftCol};
        centers = centers(:);
        centers = centers(~isnan(centers) & isfinite(centers));
        centers = round(centers);

        if isempty(centers)
            continue;
        end

        pyrCounts = extractCountsFromCenters(pyrFRs, centers, tOffsets, dtSec);
        intCounts = extractCountsFromCenters(intFRs, centers, tOffsets, dtSec);

        pyrShiftThis = cat(1, pyrShiftThis, pyrCounts);
        intShiftThis = cat(1, intShiftThis, intCounts);
    end

    pyrShiftCountsByShift{iShift} = pyrShiftThis(:);
    intShiftCountsByShift{iShift} = intShiftThis(:);

    pyrShiftMeanByShift(iShift) = mean(pyrShiftThis(:), 'omitnan');
    intShiftMeanByShift(iShift) = mean(intShiftThis(:), 'omitnan');
end

spikeStats = struct();
spikeStats.dataFile = dataFile;
spikeStats.channelsToUse = channelsToUse;
spikeStats.dtSec = dtSec;
spikeStats.nShifts = nShifts;
spikeStats.shiftColsUsed = shiftCols;

spikeStats.pyrRealCounts = pyrReal(:);
spikeStats.intRealCounts = intReal(:);

spikeStats.pyrRealMean = mean(spikeStats.pyrRealCounts, 'omitnan');
spikeStats.intRealMean = mean(spikeStats.intRealCounts, 'omitnan');

spikeStats.pyrShiftCountsByShift = pyrShiftCountsByShift;
spikeStats.intShiftCountsByShift = intShiftCountsByShift;

spikeStats.pyrShiftMeanByShift = pyrShiftMeanByShift;
spikeStats.intShiftMeanByShift = intShiftMeanByShift;

spikeStats.pyrShiftOverRealByShift = pyrShiftMeanByShift ./ spikeStats.pyrRealMean;
spikeStats.intShiftOverRealByShift = intShiftMeanByShift ./ spikeStats.intRealMean;

spikeStats.pyrShiftMeanAcrossShifts = mean(pyrShiftMeanByShift, 'omitnan');
spikeStats.intShiftMeanAcrossShifts = mean(intShiftMeanByShift, 'omitnan');

spikeStats.pyrShiftOverRealAcrossShifts = spikeStats.pyrShiftMeanAcrossShifts / spikeStats.pyrRealMean;
spikeStats.intShiftOverRealAcrossShifts = spikeStats.intShiftMeanAcrossShifts / spikeStats.intRealMean;

fprintf('\n===== spike/activity count sanity check across shifts =====\n');
fprintf('using validTransitionsNeurShiftedCell columns %d to %d\n', shiftCols(1), shiftCols(end));
fprintf('number of shifts = %d\n\n', nShifts);

fprintf('pyramidal real mean = %.6f\n', spikeStats.pyrRealMean);
fprintf('pyramidal shift mean across shifts = %.6f\n', spikeStats.pyrShiftMeanAcrossShifts);
fprintf('pyramidal shift/real across shifts = %.6f\n', spikeStats.pyrShiftOverRealAcrossShifts);
fprintf('pyramidal shift/real range = %.6f to %.6f\n\n', ...
    min(spikeStats.pyrShiftOverRealByShift, [], 'omitnan'), ...
    max(spikeStats.pyrShiftOverRealByShift, [], 'omitnan'));

fprintf('interneuron real mean = %.6f\n', spikeStats.intRealMean);
fprintf('interneuron shift mean across shifts = %.6f\n', spikeStats.intShiftMeanAcrossShifts);
fprintf('interneuron shift/real across shifts = %.6f\n', spikeStats.intShiftOverRealAcrossShifts);
fprintf('interneuron shift/real range = %.6f to %.6f\n', ...
    min(spikeStats.intShiftOverRealByShift, [], 'omitnan'), ...
    max(spikeStats.intShiftOverRealByShift, [], 'omitnan'));

figure('Color','w','Name','Shift/Real Spike Count Ratio Across Shifts');
plot(1:nShifts, spikeStats.pyrShiftOverRealByShift, 'b-o', 'LineWidth', 1.5); hold on;
plot(1:nShifts, spikeStats.intShiftOverRealByShift, 'r-o', 'LineWidth', 1.5);
yline(1, 'k--', 'LineWidth', 1.5);
xlabel('Shift Number');
ylabel('Shifted / Real Mean Count');
title('Shifted vs Real Spike/Activity Counts Across Shifts');
legend({'Pyramidal','Interneuron','Equal Activity'}, 'Location','best');
box off;

figure('Color','w','Name','Pyramidal Shift/Real Ratios Across Shifts');
histogram(spikeStats.pyrShiftOverRealByShift, 30, 'EdgeColor','none');
xline(1, 'k--', 'LineWidth', 1.5);
xlabel('Shifted / Real Mean Count');
ylabel('Number of Shifts');
title(sprintf('Pyramidal Shift/Real Ratios Across Shifts, mean = %.3f', ...
    spikeStats.pyrShiftOverRealAcrossShifts));
box off;

figure('Color','w','Name','Interneuron Shift/Real Ratios Across Shifts');
histogram(spikeStats.intShiftOverRealByShift, 30, 'EdgeColor','none');
xline(1, 'k--', 'LineWidth', 1.5);
xlabel('Shifted / Real Mean Count');
ylabel('Number of Shifts');
title(sprintf('Interneuron Shift/Real Ratios Across Shifts, mean = %.3f', ...
    spikeStats.intShiftOverRealAcrossShifts));
box off;

save(fullfile(baseDir, 'spikeCountShiftSanityCheck_AllShifts.mat'), 'spikeStats', '-v7.3');

end

function counts = extractCountsFromCenters(frMat, centers, tOffsets, dtSec)
% returns events x neurons approximate spike counts for windows centered at centers

nEvents = numel(centers);
nNeurons = size(frMat, 1);
nTime = size(frMat, 2);

counts = nan(nEvents, nNeurons);

for e = 1:nEvents
    idx = centers(e) + tOffsets;
    idx = idx(idx >= 1 & idx <= nTime);

    if isempty(idx)
        continue;
    end

    counts(e,:) = sum(frMat(:, idx), 2, 'omitnan')' * dtSec;
end

end
