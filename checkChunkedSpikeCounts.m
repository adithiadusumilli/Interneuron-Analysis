function spikeStats = checkChunkedSpikeCounts(dataFile, channelsToUse)
% checks whether shifted EMG-aligned neural windows have different spike/activity counts than unshifted windows
% assuming tAxis is in ms

if nargin < 2 || isempty(channelsToUse)
    channelsToUse = 1:4;
end

S = load(dataFile, ...
    'pyrCxWinCell', 'intCxWinCell', ...
    'pyrCxWinShiftedCell', 'intCxWinShiftedCell', ...
    'tAxis');

dtSec = median(diff(S.tAxis)) / 1000;

pyrReal = [];
intReal = [];
pyrShift = [];
intShift = [];

for ch = channelsToUse
    pyrReal = cat(1, pyrReal, sum(S.pyrCxWinCell{ch}, 3, 'omitnan') * dtSec);
    intReal = cat(1, intReal, sum(S.intCxWinCell{ch}, 3, 'omitnan') * dtSec);

    pyrShift = cat(1, pyrShift, sum(S.pyrCxWinShiftedCell{ch}, 3, 'omitnan') * dtSec);
    intShift = cat(1, intShift, sum(S.intCxWinShiftedCell{ch}, 3, 'omitnan') * dtSec);
end

spikeStats = struct();
spikeStats.dataFile = dataFile;
spikeStats.channelsToUse = channelsToUse;
spikeStats.dtSec = dtSec;

spikeStats.pyrRealCounts = pyrReal(:);
spikeStats.intRealCounts = intReal(:);
spikeStats.pyrShiftCounts = pyrShift(:);
spikeStats.intShiftCounts = intShift(:);

spikeStats.pyrRealMean = mean(spikeStats.pyrRealCounts, 'omitnan');
spikeStats.pyrShiftMean = mean(spikeStats.pyrShiftCounts, 'omitnan');
spikeStats.intRealMean = mean(spikeStats.intRealCounts, 'omitnan');
spikeStats.intShiftMean = mean(spikeStats.intShiftCounts, 'omitnan');

spikeStats.pyrShiftOverReal = spikeStats.pyrShiftMean / spikeStats.pyrRealMean;
spikeStats.intShiftOverReal = spikeStats.intShiftMean / spikeStats.intRealMean;

fprintf('\n===== spike/activity count sanity check =====\n');
fprintf('pyramidal real mean  = %.6f\n', spikeStats.pyrRealMean);
fprintf('pyramidal shift mean = %.6f\n', spikeStats.pyrShiftMean);
fprintf('pyramidal shift/real = %.6f\n\n', spikeStats.pyrShiftOverReal);

fprintf('interneuron real mean  = %.6f\n', spikeStats.intRealMean);
fprintf('interneuron shift mean = %.6f\n', spikeStats.intShiftMean);
fprintf('interneuron shift/real = %.6f\n', spikeStats.intShiftOverReal);

figure('Color','w','Name','Spike Count Check: Pyramidal');
histogram(spikeStats.pyrRealCounts, 50, 'FaceAlpha', 0.5, 'EdgeColor','none'); hold on;
histogram(spikeStats.pyrShiftCounts, 50, 'FaceAlpha', 0.5, 'EdgeColor','none');
xlabel('Approx. spikes per event-neuron window');
ylabel('Count');
title(sprintf('Pyramidal Counts: shifted/real = %.3f', spikeStats.pyrShiftOverReal));
legend({'Real','Shifted'}, 'Location','best');
box off;

figure('Color','w','Name','Spike Count Check: Interneuron');
histogram(spikeStats.intRealCounts, 50, 'FaceAlpha', 0.5, 'EdgeColor','none'); hold on;
histogram(spikeStats.intShiftCounts, 50, 'FaceAlpha', 0.5, 'EdgeColor','none');
xlabel('Approx. spikes per event-neuron window');
ylabel('Count');
title(sprintf('Interneuron Counts: shifted/real = %.3f', spikeStats.intShiftOverReal));
legend({'Real','Shifted'}, 'Location','best');
box off;

save(fullfile(fileparts(dataFile), 'spikeCountShiftSanityCheck.mat'), 'spikeStats');
end
