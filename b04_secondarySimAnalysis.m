% Simulation where the reference tissue parameters vary.

% Estimated runtime: 1400 seconds

addpath('./mfiles')
clearvars

outDir = './data/';

if ~exist(outDir,'dir')
    mkdir(outDir);
end
%%
listSigmaC = [0.02];
rng(12345)
load('./data/simMap.mat');
TRes = [15]; %30;

repF = 1000;

ktRR = 0.03:0.01:0.11;
veRR = 0.1:0.01:0.18;
kepRR = ktRR ./ veRR;

[sT sX sY] = size(simCt);
nVox = sX*sY;
CtClean = reshape(simCt,[sT sX*sY]);
CpClean = Cp;

nKtRR = length(ktRR);
nVeRR = length(veRR);
%%
trueKtRR = zeros(nKtRR,1);
trueVeRR = trueKtRR;
pkETM = zeros(nKtRR, nVox, 3, repF);
pkERRM = zeros(nVox,5);
pkCERRM = zeros(nKtRR, nVox, 5, repF);
estKtRR = zeros(nKtRR, repF);
estKtRRD = estKtRR;
estKepRRS = estKtRR;
tic
for nk=1:nKtRR
    nv = 5;
    [nk nv]
    CrrClean = ToftsKety(Cp,[ktRR(nk) ktRR(nk)./veRR(nv)],t);
    trueKtRR(nk)=ktRR(nk);
    trueVeRR(nk)=veRR(nv);
for p=1:repF 
    sigmaC = listSigmaC;
    Ct = CtClean + sigmaC*randn(size(CtClean));
    Cp = CpClean + sigmaC*randn(size(CpClean));
    Crr = CrrClean + 0.1*sigmaC*randn(size(CrrClean));
    %%
    i=1;
    dFactor=TRes(i)/initTRes;
    phaseValues = randi([0 dFactor-1],nVox,1);
    
    % Do TM, ETM
    for j=1:nVox
        curT = downsample(t, dFactor,phaseValues(j));
        curCt = downsample(Ct(:,j), dFactor,phaseValues(j));
        curCp = downsample(Cp, dFactor, phaseValues(j));
        curCrr = downsample(Crr, dFactor, phaseValues(j));
        pkETM(nk,j,:,p)=Tofts_LLSQ(curCt,curCp,curT,1);
        pkERRM(j,:)=ERRM(curCt,curCrr,curT);
    end
    
    % Get estimate for kepRR - from ERRM
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

    % Do CERRM, CLRRM, and RRIFT
    for j=1:nVox
        curT = downsample(t, dFactor,phaseValues(j));
        curCt = downsample(Ct(:,j), dFactor,phaseValues(j));
        curCrr = downsample(Crr, dFactor, phaseValues(j));
        pkCERRM(nk,j,:,p) = CERRM(curCt,curCrr,curT,estKepRR);
    end
    curT = downsample(t, dFactor,phaseValues(1));
    curCrr = downsample(Crr, dFactor, phaseValues(1));
    curCp = downsample(Cp, dFactor, phaseValues(1));
    fTail = find(curT>3,1);
    estKtRR(nk,p) = RRIFT(curCp(fTail:end),curCrr(fTail:end),curT(fTail:end),estKepRR);
    estKtRRD(nk,p) = RRIFT_diff(curCp(fTail:end),curCrr(fTail:end),curT(fTail:end),estKepRR);
    estKepRRS(nk,p) = estKepRR;
end
end
%%
trueKepRR = trueKtRR./trueVeRR;
outFile = fullfile(outDir,'simResultsTRes15-varKtRR.mat');
save(outFile,'pkETM','pkCERRM','estKtRR','estKtRRD','estKepRRS',...
    'trueKtRR','trueVeRR','trueKepRR','kepRR','ktRR','veRR',...
    'listSigmaC','TRes','t','repF');
%%
trueKtRR = zeros(nVeRR,1);
trueVeRR = trueKtRR;
pkETM = zeros(nVeRR, nVox, 3, repF);
pkERRM = zeros(nVox,5);
pkCERRM = zeros(nVeRR, nVox, 5, repF);
estKtRR = zeros(nVeRR, repF);
estKtRRD = estKtRR;
estKepRRS = estKtRR;

nk=5;
for nv=1:nVeRR
    [nk nv]
    CrrClean = ToftsKety(Cp,[ktRR(nk) ktRR(nk)./veRR(nv)],t);
    trueKtRR(nv)=ktRR(nk);
    trueVeRR(nv)=veRR(nv);
for p=1:repF 
    sigmaC = listSigmaC;
    Ct = CtClean + sigmaC*randn(size(CtClean));
    Cp = CpClean + sigmaC*randn(size(CpClean));
    Crr = CrrClean + 0.1*sigmaC*randn(size(CrrClean));
    %%
    i=1;
    dFactor=TRes(i)/initTRes;
    phaseValues = randi([0 dFactor-1],nVox,1);
    
    % Do TM, ETM
    for j=1:nVox
        curT = downsample(t, dFactor,phaseValues(j));
        curCt = downsample(Ct(:,j), dFactor,phaseValues(j));
        curCp = downsample(Cp, dFactor, phaseValues(j));
        curCrr = downsample(Crr, dFactor, phaseValues(j));
        pkETM(nv,j,:,p)=Tofts_LLSQ(curCt,curCp,curT,1);
        pkERRM(j,:)=ERRM(curCt,curCrr,curT);
    end
    
    % Get estimate for kepRR - from ERRM
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

    % Do CERRM, CLRRM, and RRIFT
    for j=1:nVox
        curT = downsample(t, dFactor,phaseValues(j));
        curCt = downsample(Ct(:,j), dFactor,phaseValues(j));
        curCrr = downsample(Crr, dFactor, phaseValues(j));
        pkCERRM(nv,j,:,p) = CERRM(curCt,curCrr,curT,estKepRR);
    end
    curT = downsample(t, dFactor,phaseValues(1));
    curCrr = downsample(Crr, dFactor, phaseValues(1));
    curCp = downsample(Cp, dFactor, phaseValues(1));
    fTail = find(curT>3,1);
    estKtRR(nv,p) = RRIFT(curCp(fTail:end),curCrr(fTail:end),curT(fTail:end),estKepRR);
    estKtRRD(nv,p) = RRIFT_diff(curCp(fTail:end),curCrr(fTail:end),curT(fTail:end),estKepRR);
    estKepRRS(nv,p) = estKepRR;
end
end
toc
%%
trueKepRR = trueKtRR./trueVeRR;
outFile = fullfile(outDir,'simResultsTRes15-varVeRR.mat');
save(outFile,'pkETM','pkCERRM','estKtRR','estKtRRD','estKepRRS',...
    'trueKtRR','trueVeRR','trueKepRR','kepRR','ktRR','veRR',...
    'listSigmaC','TRes','t','repF');
