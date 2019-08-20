%% Description
% - Loads formatted data after c01_preprocess.m
% - Fits reference region model
% - Applies reference region and input function tail (RRIFT) to estimate
%   reference tissue parameters
% - Fits Tofts model
% - Saves results to file

% Estimated runtime: ~20 seconds

%% Initialize 
clearvars
fclose('all');
addpath('./mfiles')

%% Configuration

inDir = './data/TCGA-GBM-Results/c01_preprocessed'; % Input data (from previous script)
outDir = './data/TCGA-GBM-Results/c02_postprocessed'; % Output directory

doOverwrite = true;

%% Main

if ~exist(outDir,'dir')
    mkdir(outDir)
end

matFiles = dir([inDir '/*.mat']);

% Stores R^2 of the RRIFT fit for the measured AIF tail and
% population-based AIF tail
[Rsq, RsqPop] = deal(zeros(length(matFiles),1));

tic;
for i=1:length(matFiles)
    curFile = matFiles(i).name;
    outFile = fullfile(outDir,curFile);
    if exist(outFile,'file') && ~doOverwrite
        continue
    end
    
    % Display patient ID in console as a way of tracking progress
    curFile
    
    %% Load data
    load(fullfile(inDir,curFile)); 
    % Provides: 'Ct','Cp','Crr','t','maskCt','maskCrr','maskCp' ... and more
    % refer to file b01 for details
    
    %% Basic pre-processing
    
    % Set negative concentrations to zero
    Crr(Crr<0)=0;
    Ct(Ct<0)=0;
    Cp(Cp<0)=0;
    
    % Identify voxels with negligible enhancement.
    % In this case, 'negligible' is when the maximum concentration is below 0.01 mM.
    enhancementMask = max(Ct) > 0.01;
    % Keep track of how many voxels there were
    numGoodVox = sum(enhancementMask(:));
    numVox = sum(maskCt(:));
    % Remove voxels with negligible enhancement
    Ct = Ct(:,enhancementMask);
    maskCt(maskCt) = enhancementMask;
    
    %% Fit Reference Region Model
    [pkCE, ~, estKepRR, pkERRM] = CERRM(Ct,Crr,t);
    % `estKepRR` will be used by RRIFT in next block
    %% Apply RRIFT
    % Tail begins at ~33rd frame since it is ~3 minutes into acquisition
    % or ~2 minutes after bolus arrival since bolus arrival is at ~1 minute
    % Alternative to hard-cording the frame namer: find(t>3, 1, 'first')
    fTail = 33;
    tailCrr = Crr(fTail:end);
    tailT = t(fTail:end);
    
    % RRIFT - with measured AIF tail
    tailCp = Cp(fTail:end);
    [estKtRR, num, denum] = RRIFT(tailCp,tailCrr,tailT,estKepRR);
    estVeRR = estKtRR./estKepRR;
    Rsq(i) = corr2(num, denum)^2;
    
    % RRIFT - with assumed population-based tail
    CpPopAvg = GeorgiouAif(t,t(7));
    tailCpPop = CpPopAvg(fTail:end);
    [estKtRRPop, numPop, denumPop] = RRIFT(tailCpPop,tailCrr,tailT,estKepRR);
    estVeRRPop = estKtRRPop./estKepRR;
    RsqPop(i) = corr2(numPop, denumPop)^2;
    
    % Alternative RRIFT - using differentiation instead of integration [EXTRA]
    [estKtRRdiff] = RRIFT_diff(tailCp,tailCrr,tailT,estKepRR);
    estVeRRdiff = estKtRRdiff./estKepRR;

    %% Fit Tofts Model
    ETM = struct;
    ETM.tumour = Tofts_LLSQ(Ct,Cp,t,1); % Extended Tofts model fit on tumour voxels
    ETM.muscle = Tofts_LLSQ(Crr,Cp,t,1); % Extended Tofts model fit on muscle region 
    
    %% Estimate fitting uncertainty in muscle (for error-bars)
    % For RRIFT - KepRR
    rawKepRR = pkERRM(:,5);
    goodVals = pkERRM(:,1)>0 & pkERRM(:,2)>0 & pkERRM(:,3)>0 & pkERRM(:,4)>0 & pkERRM(:,5)>0 & imag(pkERRM(:,5))==0;
    rawKepRR = rawKepRR(goodVals);
    qtRange = quantile(rawKepRR,[.25 .75]);
    kepMask = rawKepRR>qtRange(1) & rawKepRR<qtRange(2);
    ciKepRR = peak2peak(ConfInterval(rawKepRR(kepMask)'));
    % For RRIFT - KtransRR
    [ahat,r,J] = nlinfit(denum,num,@(a,x)(a*x),[0]);
    ciKtRR = peak2peak(nlparci(ahat,r,'Jacobian',J));
    % For RRIFT - veRR = KtransRR/veRR
    ciVeRR = sqrt(estVeRR^2 * ( (ciKtRR/estKtRR)^2 + (ciKepRR/estKepRR)^2))/2;
    % For Tofts Model
    [pkE,r,J] = nlinfit(t,Crr,@(a,x)ToftsKety(Cp, a, x),[0.1 0.1 0.1]);
    cis = nlparci(pkE,r,'Jacobian',J);
    ciKepRRT = peak2peak(cis(2,:));
    ciKtRRT = peak2peak(cis(1,:));
    estVeRRT = pkE(1)./pkE(2);
    ciVeRRT = sqrt(estVeRRT^2 * ( (ciKtRRT/pkE(1))^2 + (ciKepRRT/pkE(2))^2))/2;
    %% Make the maps
    [sX, sY, sZ] = size(maskCt);
    [mapKt, mapKep, mapVe, mapVp, mapKtR, mapKepR, mapVeR, mapVpR] = ...
        deal(zeros(sX,sY,sZ));
    
    % Maps from extended Tofts model
    mapKt(maskCt)=ETM.tumour(:,1);
    mapKep(maskCt)=ETM.tumour(:,2);
    mapVe(maskCt)=mapKt(maskCt)./mapKep(maskCt);
    mapVp(maskCt)=ETM.tumour(:,3);

    % Maps from RRIFT using the measured AIF tail
    mapKtR(maskCt)=estKtRR * pkCE(:,1);
    mapKepR(maskCt)=pkCE(:,3);
    mapVeR(maskCt)=estVeRR * pkCE(:,2);
    mapVpR(maskCt)=estKtRR * pkCE(:,4);
    
    %% Save results
    
    save(outFile,...
        'mapKt','mapKep','mapVe','mapVp',...
        'mapKtR','mapKepR','mapVeR','mapVpR',...
        'Crr','Cp','t','maskCt','maskCrr','ETM','estKepRR',...
        'estKtRR','estVeRR','estKtRRPop','estVeRRPop','estKtRRdiff','estVeRRdiff',...
        'Rsq','RsqPop','num','denum','numVox','numGoodVox',...
        'cnr','snr','sigmaCt','T1Cp','T1Crr','T1Ct',...
        'ciKtRR','ciKepRR','ciVeRR','ciKtRRT','ciKepRRT','ciVeRRT');
end
toc