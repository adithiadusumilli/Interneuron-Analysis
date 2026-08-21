function combineChunkedPopavgBehaviorOutputsAllSessions()
% Combines finalized behavior-specific chunked trial-averaged popavg Quest outputs.
clc;

baseDirs = {'/home/asa7288/Transfer/D026','/home/asa7288/Transfer/D020', ...
    '/home/asa7288/Transfer/D024','/home/asa7288/Transfer/D043', ...
    '/home/asa7288/Transfer/D050','/home/asa7288/Transfer/D054'};
mouseIDs = {'D026','D020','D024','D043','D050','D054'};
labelsToUse = 'classifier'; labelTag = lower(labelsToUse);
behaviors = 1:10; nPermExpected = 300; nShiftExpected = 100;
outFile = '/home/asa7288/concatCrossCorrPerCanonicalBehavior_classifier_chunked_trialavg_popavg_ALL_SESSIONS.mat';

allSessions = struct('baseDirs',{baseDirs},'mouseIDs',{mouseIDs}, ...
    'labelsToUse',labelsToUse,'behaviors',behaviors, ...
    'nPermExpected',nPermExpected,'nShiftExpected',nShiftExpected);
allSessions.permBlocks = struct('H0',1:100,'H50',101:200,'Hneg50',201:300);

for s = 1:numel(baseDirs)
    baseDir = baseDirs{s}; outDir = fullfile(baseDir,'quest_runs');
    sess = struct('baseDir',baseDir,'mouseID',mouseIDs{s},'sessionTag',mouseIDs{s}, ...
        'labelsToUse',labelsToUse,'behaviors',[]);
    fprintf('\n================ %s ================\n',mouseIDs{s});

    for bi = 1:numel(behaviors)
        beh = behaviors(bi);
        B = struct();
        B.beh = beh;

        %% REAL
        rf = fullfile(outDir,sprintf( ...
            'concatCrossCorrPerCanonicalBehavior_%s_unperm_beh%02d.mat',labelTag,beh));
        B.real_file_found = isfile(rf);
        B.real_xc=[]; B.real_peakLagSec=NaN; B.real_peakCorr=NaN;
        B.real_nTrialsThisBehavior=NaN; B.real_nTrialsUsed=NaN;
        B.real_pyrAvgTrace=[]; B.real_intAvgTrace=[];
        B.lags=[]; B.binSize=NaN; B.chunkHalf=NaN;
        B.doBaselineNorm=[]; B.channelsToUseThisSess=[];
        if B.real_file_found
            R=load(rf);
            B.real_xc=g(R,'xc',[]); B.real_peakLagSec=g(R,'peakLagSec',NaN);
            B.real_peakCorr=g(R,'peakCorr',NaN);
            B.real_nTrialsThisBehavior=g(R,'nTrialsThisBehavior',NaN);
            B.real_nTrialsUsed=g(R,'nTrialsUsed',NaN);
            B.real_pyrAvgTrace=g(R,'pyrAvgTrace',[]);
            B.real_intAvgTrace=g(R,'intAvgTrace',[]);
            B.lags=g(R,'lags',[]); B.binSize=g(R,'binSize',NaN);
            B.chunkHalf=g(R,'chunkHalf',NaN);
            B.doBaselineNorm=g(R,'doBaselineNorm',[]);
            B.channelsToUseThisSess=g(R,'channelsToUseThisSess',[]);
        end

        %% PERMUTATIONS -- fixed row = permutation index
        B.perm_peakLagSec=nan(nPermExpected,1);
        B.perm_peakCorr=nan(nPermExpected,1);
        B.perm_accepted=false(nPermExpected,1);
        B.perm_tryCount=nan(nPermExpected,1);
        B.perm_nTrialsThisBehavior=nan(nPermExpected,1);
        B.perm_nTrialsUsed=nan(nPermExpected,1);
        B.perm_runtimeSec=nan(nPermExpected,1);
        B.perm_savedShiftCorrUpper95=nan(nPermExpected,1);
        B.perm_file_found=false(nPermExpected,1);
        B.perm_xc=[];

        for p=1:nPermExpected
            pf=fullfile(outDir,sprintf( ...
                'concatCrossCorrPerCanonicalBehavior_%s_perm_%03d_beh%02d.mat', ...
                labelTag,p,beh));
            if ~isfile(pf), continue; end
            P=load(pf); B.perm_file_found(p)=true;
            B.perm_peakLagSec(p)=g(P,'peakLagSec',NaN);
            B.perm_peakCorr(p)=g(P,'peakCorr',NaN);
            B.perm_accepted(p)=logical(g(P,'acceptedPerm',false));
            B.perm_tryCount(p)=g(P,'tryCount',NaN);
            B.perm_nTrialsThisBehavior(p)=g(P,'nTrialsThisBehavior',NaN);
            B.perm_nTrialsUsed(p)=g(P,'nTrialsUsed',NaN);
            B.perm_runtimeSec(p)=g(P,'runtimeSec',NaN);
            B.perm_savedShiftCorrUpper95(p)=g(P,'shiftCorrUpper95',NaN);
            if isempty(B.lags), B.lags=g(P,'lags',[]); end
            if isnan(B.binSize), B.binSize=g(P,'binSize',NaN); end
            if isnan(B.chunkHalf), B.chunkHalf=g(P,'chunkHalf',NaN); end
            if isempty(B.doBaselineNorm), B.doBaselineNorm=g(P,'doBaselineNorm',[]); end
            if isempty(B.channelsToUseThisSess)
                B.channelsToUseThisSess=g(P,'channelsToUseThisSess',[]);
            end
            if isfield(P,'xc') && ~isempty(P.xc)
                if isempty(B.perm_xc)
                    B.perm_xc=nan(nPermExpected,numel(P.xc));
                end
                if size(B.perm_xc,2)==numel(P.xc), B.perm_xc(p,:)=P.xc(:)'; end
            end
        end
        B.perm_inds_found=find(B.perm_file_found);
        B.perm_files_expected=nPermExpected;
        B.perm_files_found=nnz(B.perm_file_found);
        B.perm_files_accepted=nnz(B.perm_file_found & B.perm_accepted);

        %% SHIFTS -- fixed row = shift index
        B.shift_xcZeroLag=nan(nShiftExpected,1);
        B.shift_nTrialsThisBehavior=nan(nShiftExpected,1);
        B.shift_nTrialsUsed=nan(nShiftExpected,1);
        B.shift_usedShiftCol=nan(nShiftExpected,1);
        B.shift_usedShiftAmtMs=nan(nShiftExpected,1);
        B.shift_usedFallbackShift=false(nShiftExpected,1);
        B.shift_retryCount=nan(nShiftExpected,1);
        B.shift_file_found=false(nShiftExpected,1);

        for sh=1:nShiftExpected
            sf=fullfile(outDir,sprintf( ...
                'concatCrossCorrPerCanonicalBehavior_%s_shift_%03d_zerolag_beh%02d.mat', ...
                labelTag,sh,beh));
            if ~isfile(sf), continue; end
            D=load(sf); B.shift_file_found(sh)=true;
            B.shift_xcZeroLag(sh)=g(D,'xcZeroLag',NaN);
            B.shift_nTrialsThisBehavior(sh)=g(D,'nTrialsThisBehavior',NaN);
            B.shift_nTrialsUsed(sh)=g(D,'nTrialsUsed',NaN);
            B.shift_usedShiftCol(sh)=g(D,'usedShiftCol',NaN);
            B.shift_usedShiftAmtMs(sh)=g(D,'usedShiftAmtMs',NaN);
            B.shift_usedFallbackShift(sh)=logical(g(D,'usedFallbackShift',false));
            B.shift_retryCount(sh)=g(D,'retryCount',NaN);
        end
        B.shift_inds_found=find(B.shift_file_found);
        B.shift_files_expected=nShiftExpected;
        B.shift_files_found=nnz(B.shift_file_found);
        vals=B.shift_xcZeroLag(isfinite(B.shift_xcZeroLag));
        B.nFiniteShiftControls=numel(vals);
        if isempty(vals), B.shiftCorrUpper95=NaN;
        else, B.shiftCorrUpper95=prctile(vals,95); end

        % Sanity check against thresholds saved inside permutation files.
        saved=B.perm_savedShiftCorrUpper95(isfinite(B.perm_savedShiftCorrUpper95));
        if ~isempty(saved) && isfinite(B.shiftCorrUpper95)
            d=max(abs(saved-B.shiftCorrUpper95));
            if d>1e-12
                warning('%s beh %d threshold mismatch; max diff %.12g',mouseIDs{s},beh,d);
            end
        end

        fprintf('beh %2d | real %d | perms %3d/300 (%3d accepted) | shifts %3d/100 | thresh %.6f\n', ...
            beh,B.real_file_found,B.perm_files_found,B.perm_files_accepted, ...
            B.shift_files_found,B.shiftCorrUpper95);
        sess.behaviors(bi)=B;
    end
    allSessions.sessions(s)=sess;
end

%% SUMMARY TABLE
n=numel(baseDirs)*numel(behaviors);
Animal=strings(n,1); Behavior=nan(n,1); RealFound=false(n,1);
RealPeakLag_s=nan(n,1); RealPeakCorr=nan(n,1);
PermFilesFound=nan(n,1); PermAccepted=nan(n,1);
ShiftFilesFound=nan(n,1); FiniteShiftControls=nan(n,1);
ShiftCorrUpper95=nan(n,1); r=0;
for s=1:numel(baseDirs)
    for bi=1:numel(behaviors)
        r=r+1; B=allSessions.sessions(s).behaviors(bi);
        Animal(r)=string(mouseIDs{s}); Behavior(r)=B.beh;
        RealFound(r)=B.real_file_found; RealPeakLag_s(r)=B.real_peakLagSec;
        RealPeakCorr(r)=B.real_peakCorr; PermFilesFound(r)=B.perm_files_found;
        PermAccepted(r)=B.perm_files_accepted; ShiftFilesFound(r)=B.shift_files_found;
        FiniteShiftControls(r)=B.nFiniteShiftControls;
        ShiftCorrUpper95(r)=B.shiftCorrUpper95;
    end
end
summaryTable=table(Animal,Behavior,RealFound,RealPeakLag_s,RealPeakCorr, ...
    PermFilesFound,PermAccepted,ShiftFilesFound,FiniteShiftControls,ShiftCorrUpper95);
allSessions.summaryTable=summaryTable;

save(outFile,'allSessions','summaryTable','-v7.3');
fprintf('\nSaved combined master file:\n%s\n',outFile);
disp(summaryTable);
end

function v=g(S,name,defaultVal)
if isfield(S,name), v=S.(name); else, v=defaultVal; end
end
