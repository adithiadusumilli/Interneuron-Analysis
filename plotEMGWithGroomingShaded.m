function plotEMGWithGroomingShaded(baseDirs)

% plots raw EMG for each animal with grooming epochs shaded
% this makes it easy to visually compare EMG during grooming

behaviorToPlot = 4; % grooming
nAnimals = numel(baseDirs);
nCh = 4;

figure('Color','w')
tiledlayout(nAnimals,1,'Padding','compact','TileSpacing','compact')

for iDir = 1:nAnimals
    
    baseDir = baseDirs{iDir};
    
    % extract animal name
    tok = regexp(baseDir,'D\d+','match','once');
    if isempty(tok)
        animalID = sprintf('animal%d',iDir);
    else
        animalID = tok;
    end
    
    % load data
    E = load(fullfile(baseDir,'EMG1ms.mat'),'downsampEMG');
    U = load(fullfile(baseDir,'UMAP.mat'), ...
        'origDownsampEMGInd','classifierLabels','classifierBehvs');
    
    emg = E.downsampEMG(1,:); % just channel 1 for quick diagnostic
    nTime = numel(emg);
    
    % map classifier labels back to full EMG timeline
    labels1k = nan(1,nTime);
    
    origInd = U.origDownsampEMGInd(:);
    lab = U.classifierLabels(:);
    
    labels1k(origInd) = lab;
    
    % grooming mask
    mask = labels1k == behaviorToPlot;
    
    nexttile
    hold on
    
    t = (0:nTime-1)/1000; % seconds
    
    plot(t,emg,'k')
    
    % shade grooming
    d = diff([0 mask 0]);
    s = find(d==1);
    e = find(d==-1)-1;
    
    for k = 1:numel(s)
        xs = t([s(k) e(k)]);
        patch([xs(1) xs(2) xs(2) xs(1)], ...
              [min(emg) min(emg) max(emg) max(emg)], ...
              [1 0.8 0.8], ...
              'EdgeColor','none','FaceAlpha',0.3)
    end
    
    plot(t,emg,'k')
    
    title(sprintf('%s EMG (grooming shaded)',animalID))
    xlabel('time (s)')
    ylabel('EMG')
    
end
end
