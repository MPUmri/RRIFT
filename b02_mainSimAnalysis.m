% Improved version of the simMap downsample which also adds noise and
% replicates each compartment

% Estimated runtime: 1700 seconds

addpath('./mfiles')
clearvars

outDir = './data/';

if ~exist(outDir,'dir')
    mkdir(outDir);
end
%%
rng(12345) % Fix the random seed for reproducibility
load('./data/simMap.mat');

TRes = [5,10,15,30]; % range of temporal resolutions, in seconds
listSigmaC = 0:0.01:0.05; % range of noise levels, in mM
repF = 1000; % Number of replications for each scenario

% Reference tissue properties
ktRR = 0.07;
kepRR = 0.5;
veRR = ktRR/kepRR;
Crr = ToftsKety(Cp,[ktRR,kepRR],t);

[sT, sX, sY] = size(simCt);
nVox = sX*sY;
CtClean = reshape(simCt,[sT sX*sY]);
CpClean = Cp;
CrrClean = Crr;
%%
params = struct;
params.ETM = zeros(nVox,3,repF,length(TRes),length(listSigmaC));
params.CERRM = zeros(nVox,5,repF,length(TRes),length(listSigmaC));

pkERRM = zeros(nVox,5);

estKtRR = zeros(repF,length(TRes),length(listSigmaC));
estKtRRD = estKtRR;
estKepRRS = estKtRR;

hWait = waitbar(0,'Running simulation...');

tic
for p=1:repF
    waitbar(p/repF, hWait, sprintf('Running simulation [%d/%d]...', p, repF))
for q=1:length(listSigmaC)
    sigmaC = listSigmaC(q);

    Ct = CtClean + sigmaC*randn(size(CtClean));
    Cp = CpClean + sigmaC*randn(size(CpClean));
    Crr = CrrClean + 0.1*sigmaC*randn(size(CrrClean));
for i=1:length(TRes)
    %%
    dFactor=TRes(i)/initTRes;
    phaseValues = randi([0 dFactor-1],nVox,1);
    
    % Do TM, ETM
    for j=1:nVox
        curT = downsample(t, dFactor,phaseValues(j));
        curCt = downsample(Ct(:,j), dFactor,phaseValues(j));
        curCp = downsample(Cp, dFactor, phaseValues(j));
        curCrr = downsample(Crr, dFactor, phaseValues(j));
        params.ETM(j,:,p,i,q)=Tofts_LLSQ(curCt,curCp,curT,1);
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
        params.CERRM(j,:,p,i,q) = CERRM(curCt,curCrr,curT,estKepRR);
    end
    curT = downsample(t, dFactor,phaseValues(1));
    curCrr = downsample(Crr, dFactor, phaseValues(1));
    curCp = downsample(Cp, dFactor, phaseValues(1));
    fTail = find(curT>3,1);
    estKtRR(p,i,q) = RRIFT(curCp(fTail:end), curCrr(fTail:end), curT(fTail:end), estKepRR);
    estKtRRD(p,i,q) = RRIFT_diff(curCp(fTail:end), curCrr(fTail:end), curT(fTail:end), estKepRR);
    estKepRRS(p,i,q) = estKepRR;
    %% Look at how number of tail frames affects results
    if TRes(i) == 5
        zind=1;
        for z = fTail:length(curT)-1
            tailT_5(p,q,zind) = curT(z);
            estKtRRS_5(p,q,zind) = RRIFT(curCp(z:end), curCrr(z:end), curT(z:end), estKepRR);
            zind = zind+1;
        end
    elseif TRes(i) == 15
        zind=1;
        for z = fTail:length(curT)-1
            tailT_15(p,q,zind) = curT(z);
            estKtRRS_15(p,q,zind) = RRIFT(curCp(z:end), curCrr(z:end), curT(z:end), estKepRR);
            zind = zind+1;
        end
    elseif TRes(i) == 30
        zind=1;
        for z = fTail:length(curT)-1
            tailT_30(p,q,zind) = curT(z);
            estKtRRS_30(p,q,zind) = RRIFT(curCp(z:end), curCrr(z:end), curT(z:end), estKepRR);
            zind = zind+1;
        end
    end
end
end
end
close(hWait)
%%
outFile = fullfile(outDir,'simResults.mat');
save(outFile,'params','estKtRR','estKtRRD','estKepRRS','kepRR','ktRR','veRR',...
    'listSigmaC','TRes','t','repF',...
    'tailT_5','tailT_15','tailT_30','estKtRRS_5','estKtRRS_15','estKtRRS_30');
toc