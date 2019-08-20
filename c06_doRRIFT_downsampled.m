% This script applies RRIFT and the ETM on downsampled in-vivo data.

% Estimated runtime: 600 seconds

clearvars
fclose('all')
addpath('./mfiles')

inDir = './data/TCGA-GBM-Results/c01_preprocessed';

outMapDir = './data/TCGA-GBM-Results/c06_downsampled';

if ~exist(outMapDir)
    mkdir(outMapDir)
end

matFiles = dir([inDir '/*.mat']);

tic
for q=1:length(matFiles)
    curFile = matFiles(q).name;
    outFile = fullfile(outMapDir,curFile);

    load(fullfile(inDir,curFile));
    % 'Ct','Cp','Crr','t','maskCt','maskCrr','maskCp'
    Crr(Crr<0)=0;
    Ct(Ct<0)=0;
    CpPop = GeorgiouAif(t,t(7));
    Cp(Cp<0)=0;
    
    qtMask = max(Ct) > 0.01;
    numGoodVox = sum(qtMask(:));
    numVox = sum(maskCt(:));
    Ct = Ct(:,qtMask);
    maskCt(maskCt) = qtMask;
    %%
    numVox = sum(maskCt(:));
    dFactors = 1:10;
    pkETM = zeros(numVox,3,length(dFactors));
    pkCE = zeros(numVox,5,length(dFactors));
for i=1:length(dFactors)
    dFactor = dFactors(i);
    phaseValues = randi([0 dFactor-1],numVox,1);
    
    pkERRM = zeros(numVox,5);
    for j=1:numVox
        curT = downsample(t,dFactor, phaseValues(j));
        curCt = downsample(Ct(:,j), dFactor, phaseValues(j));
        curCp = downsample(Cp, dFactor, phaseValues(j));
        curCrr = downsample(Crr, dFactor, phaseValues(j));
        pkETM(j,:,i)=Tofts_LLSQ(curCt,curCp,curT,1);
        pkERRM(j,:)=ERRM(curCt,curCrr,curT);
    end

    rawKepRR = pkERRM(:,5);
    if std(rawKepRR)<1e-3
       % If the estimated kepRR from ERRM is closely grouped, then use median
       % This situation is very unlikely to happen in clinical data
       % The practical purpose of this is that when simulating
       % noiseless data, the interquartile mean can't be used because
       % there is no fluctuation in the estimated kepRR from ERRM
       estKepRR = nanmedian(rawKepRR);
    else
       % Find voxels where all estimates are real and positive
       goodVals = pkERRM(:,1)>0 & pkERRM(:,2)>0 & pkERRM(:,3)>0 & pkERRM(:,4)>0 & pkERRM(:,5)>0 & imag(pkERRM(:,5))==0;
       estKepRR = iqrMean(rawKepRR(goodVals));
    end

    for j=1:numVox
        curT = downsample(t,dFactor, phaseValues(j));
        curCt = downsample(Ct(:,j), dFactor, phaseValues(j));
        curCp = downsample(Cp, dFactor, phaseValues(j));
        curCrr = downsample(Crr, dFactor, phaseValues(j));
        pkCE(j,:,i) = CERRM(curCt,curCrr,curT,estKepRR);
    end
    curT = downsample(t, dFactor,phaseValues(1));
    curCrr = downsample(Crr, dFactor, phaseValues(1));
    curCp = downsample(Cp, dFactor, phaseValues(1));
    curCpPop = downsample(CpPop, dFactor, phaseValues(1));
    fTail = find(curT>3,1);
    estKtRRs(i) = RRIFT(curCp(fTail:end),curCrr(fTail:end),curT(fTail:end),estKepRR);
    estKtRRsPop(i) = RRIFT(curCpPop(fTail:end),curCrr(fTail:end),curT(fTail:end),estKepRR);
    estKepRRs(i) = estKepRR;
    TRes(i) = (curT(2)-curT(1))*60;
    
    pkCE(:,1,i) = pkCE(:,1,i)*estKtRRs(i);
    pkCE(:,2,i) = pkCE(:,2,i)*(estKtRRs(i)/estKepRRs(i));
    pkCE(:,4,i) = pkCE(:,4,i)*estKtRRs(i);
end
save(fullfile(outMapDir,curFile),'dFactors','TRes','estKtRRs','estKtRRsPop','estKepRRs',...
    'maskCt','pkETM','pkCE','Crr','Cp','t');
end
toc
return
%%
errKep = PercentError(rrKep,etKep);
errKt = PercentError(rrKt,etKt);
errVe=PercentError(rrVe,etVe);
%%
plot(errKep)
%%
plot(errKt)
%%
plot(errVe)
%%
plot(1:27,errKep,1:27,errKt)
%%
gbmErrKep = errKep(1:17);
errMask = abs(gbmErrKep)<40;
goodErr = gbmErrKep(errMask);
badErr = gbmErrKep(~errMask);
scatter(1:length(goodErr), goodErr, 'filled')
hold on
scatter(1:length(badErr), badErr, 'filled')
xlim([0 12])

%%
x = abs(PercentError(rrKt,etKt)) < 40 & abs(PercentError(rrKep,etKep)) < 40;
sum(x)
%%
x = errMask;
for i=1:length(x)
    if ~x(i)
        continue
    end
    i
    h=figure;
    plot(tList{i},cpList{i})
    hold on
    plot(tList{i},estCpList{i})
    waitfor(h)
end

%%
x = errMask;
for i=1:length(x)
    if ~x(i)
        continue
    end
    i
    h=figure;
    plot(tList{i},crrList{i})
    hold on
    plot(tList{i},emmCrrList{i})
    waitfor(h)
end
%%
gVal = [2,3,4,5,6,14,16];
gVal2 = [2,3,4,5,6,14,16,7,8,13,15];
gVal=[2,3,4,5,6,7,8,10,13,14,16];

curVal = gVal;
for i=1:length(curVal)
    curI = curVal(i)
    estCrr = ToftsKety(cpList{curI},[ttKt(curI) ttKep(curI) etVp(curI)],tList{curI},0);
    h=figure;
    scatter(tList{curI},crrList{curI}); hold on;
    plot(tList{curI},estCrr)
    title(curI)
end

%%
gVal=[2,3,4,5,6,7,10,13,14,16];
scatter(ttVe(gVal),rrVe(gVal))