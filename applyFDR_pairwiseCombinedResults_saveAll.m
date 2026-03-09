function applyFDR_pairwiseCombinedResults_saveAll(saveFile, realField, nullField, lagField, xcField, alpha, flag, tailType)

% applies david's FDRcutoff.m to combined pairwise results
% saves ALL original matrices again so plotting scripts can reuse the file

if nargin < 6 || isempty(alpha)
    alpha = 0.05;
end

if nargin < 7 || isempty(flag)
    flag = true;
end

if nargin < 8 || isempty(tailType)
    tailType = "right";
end

tailType = lower(string(tailType));

S = load(saveFile,'allSessions');
allSessions = S.allSessions;

nSess = numel(allSessions.sessions);

FDRresults = struct();
FDRresults.sessions = struct([]);

for s = 1:nSess

    sess = allSessions.sessions(s);

    animalID = regexp(sess.baseDir,'D\d+','match','once');

    fprintf('\nprocessing %s\n',animalID);

    realMat = sess.(realField);
    nullMat = sess.(nullField);

    if isempty(realMat) || isempty(nullMat)
        warning('missing data, skipping session');
        continue
    end

    n = size(realMat,1);

    upperMask = triu(true(n),1);
    keepLinear = find(upperMask);

    realVals = realMat(keepLinear);

    nNull = size(nullMat,3);

    nullVals = nan(numel(keepLinear),nNull);

    for k = 1:nNull
        tmp = nullMat(:,:,k);
        nullVals(:,k) = tmp(keepLinear);
    end

    valid = ~isnan(realVals) & all(~isnan(nullVals),2);

    keepLinear = keepLinear(valid);
    realVals = realVals(valid);
    nullVals = nullVals(valid,:);

    pVals = nan(size(realVals));

    for i = 1:numel(realVals)

        switch tailType

            case "right"
                pVals(i) = (sum(nullVals(i,:) >= realVals(i)) + 1) / (nNull + 1);

            case "left"
                pVals(i) = (sum(nullVals(i,:) <= realVals(i)) + 1) / (nNull + 1);

            case "two"
                mu = mean(nullVals(i,:),'omitnan');
                realDev = abs(realVals(i)-mu);
                nullDev = abs(nullVals(i,:)-mu);
                pVals(i) = (sum(nullDev >= realDev) + 1) / (nNull + 1);

        end
    end

    fdrCut = FDRcutoff(pVals,alpha,flag);

    sigVec = pVals <= fdrCut;

    pMat = nan(n);
    sigMask = false(n);

    pMat(keepLinear) = pVals;
    sigMask(keepLinear) = sigVec;

    pMat = pMat + pMat.';
    sigMask = sigMask | sigMask.';

    %% --- save results ---

    out = struct();

    out.baseDir = sess.baseDir;
    out.animalID = animalID;

    out.fdrCutoff = fdrCut;
    out.nTests = numel(pVals);
    out.nSignificant = nnz(triu(sigMask,1));

    out.pMat = pMat;
    out.sigMaskFDR = sigMask;

    % original matrices saved again
    out.(realField) = realMat;
    out.(nullField) = nullMat;

    if isfield(sess,lagField)
        out.(lagField) = sess.(lagField);
    end

    if isfield(sess,xcField)
        out.(xcField) = sess.(xcField);
    end

    % significant-only versions
    sigReal = realMat;
    sigReal(~sigMask) = NaN;
    out.(['sig_' realField]) = sigReal;

    if isfield(sess,lagField)
        sigLag = sess.(lagField);
        sigLag(~sigMask) = NaN;
        out.(['sig_' lagField]) = sigLag;
    end

    if isfield(sess,xcField)
        sigXC = sess.(xcField);
        sigXC(~sigMask) = NaN;
        out.(['sig_' xcField]) = sigXC;
    end

    FDRresults.sessions(s) = out;

    fprintf('%s | tests=%d | significant=%d\n',animalID,out.nTests,out.nSignificant)

end

%% save file

[pathstr,name] = fileparts(saveFile);

outFile = fullfile(pathstr,[name '_FDR.mat']);

save(outFile,'FDRresults','-v7.3')

fprintf('\nsaved:\n%s\n',outFile)

end
