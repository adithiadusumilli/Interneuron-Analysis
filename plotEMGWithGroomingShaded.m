function plotEMGWithGroomingShaded(baseDirs)

% plots raw EMG for each animal with grooming epochs shaded
% this makes it easy to visually compare EMG during grooming

behaviorToPlot = 4; % grooming
nAnimals = numel(baseDirs);

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

    % make lengths match
    n = min(numel(origInd), numel(lab));
    origInd = origInd(1:n);
    lab = lab(1:n);

    % keep only valid indices
    ok = origInd >= 1 & origInd <= nTime & ~isnan(lab);
    labels1k(origInd(ok)) = lab(ok);
    
    % grooming mask
    mask = labels1k == behaviorToPlot;
    
    nexttile
    hold on
    
    t = (0:nTime-1)/1000; % seconds
    
    % first plot emg lightly so patches do not hide it
    plot(t,emg,'k')
    
    % shade grooming
    d = diff([0 mask 0]);
    s = find(d==1);
    e = find(d==-1)-1;
    
    yMin = min(emg);
    yMax = max(emg);

    for k = 1:numel(s)
        xs = t([s(k) e(k)]);
        patch([xs(1) xs(2) xs(2) xs(1)], ...
              [yMin yMin yMax yMax], ...
              [1 0.8 0.8], ...
              'EdgeColor','none','FaceAlpha',0.3)
    end
    
    % redraw emg on top
    plot(t,emg,'k')
    
    title(sprintf('%s EMG channel 1 (grooming shaded)',animalID))
    xlabel('time (s)')
    ylabel('EMG')
    
    hold off
end
end
