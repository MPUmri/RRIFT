%% Description
% - Loads fitting results from c02_doRRIFT.m
% - Makes figures for analyzing results

%% Initialize
clearvars
fclose('all');
addpath('./mfiles')

%% Configuration

inDir = './data/TCGA-GBM-Results/c02_postprocessed'; % Input directory (output from c02_doRRIFT.m)

%% Main

matFiles = dir([inDir '/*.mat']);

% Initialize arrays (the variable names will be cleaned up later)
[
    arrKtETM, arrKepETM, arrVeETM, arrVpETM,... % Tumour parameters from ETM
    arrKtRRIFT, arrKepRRIFT, arrVeRRIFT, arrVpRRIFT,... % Tumour parameters from RRIFT-measuredTail
    arrKtRRIFTPop, arrVeRRIFTPop, arrVpRRIFTPop,... % Tumour parameters from RRIFT-assumedTail
    arrKtRRM, arrVeRRM, arrVpRRM,... % Tumour parameters from reference region model (i.e. relative parameters)
    crrList, cpList, tList... % Reference tissue & AIF curve for each patient
] = deal([]); % This will waste some RAM, but it makes the code simpler
% Pre-allocate more arrays
[
    arrKtETM_m, arrKepETM_m, arrVeETM_m, arrVpETM_m,... % Muscle parameters from ETM
    arrKepRRIFT_m, arrVeRRIFT_m, arrKtRRIFT_m,... % Muscle parameters from RRIFT-measuredTail
    arrKtRRIFTPop_m, arrVeRRIFTPop_m,... % Muscle parameters from RRIFT-assumedTail
    arrKtRRIFTD_m, arrVeRRIFTD_m,... % Muscle parameters from differential form of RRIFT
    cccKt, cccKep, cccVe, cccVp,...  % Concordance correlation coefficient for each study
    t1CpList, t1CrrList, t1CtList,... % T1 for the AIF, muscle, and tumour
    nVoxGood, nVox, nVoxRR, snrList, cnrList, sigmaList,... % miscellanea
    arrKtRRci, arrKepRRci, arrVeRRci, arrKtRRTci, arrKepRRTci, arrVeRRTci... % uncertainty of muscle estimates
] = deal(zeros(length(matFiles),1));

for i=1:length(matFiles)
    curFile = matFiles(i).name;
    
    load(fullfile(inDir,curFile));
    % Loads (refer to b02_doRRIFT.m):
        % mapKt, mapKep, mapVe, mapVp, ...
        % mapKtR, mapKepR, mapVeR, mapVpR, ...
        % Crr, Cp, t, maskCt, maskCrr, ETM, estKepRR, ...
        % estKtRR, estVeRR, estKtRRPop, estVeRRPop, estKtRRdiff, estVeRRdiff, ...
        % Rsq, RsqPop, num, denum, numVox, numGoodVox, ...
        % cnr, snr, sigmaCt, T1Cp, T1Crr, T1Ct, ...
        % ciKtRR, ciKepRR, ciVeRR, ciKtRRT, ciKepRRT, ciVeRRT
    %% Concatenate all the fitted parameters for tumour
    % From extended Tofts model
    arrKtETM = [arrKtETM; mapKt(maskCt)];
    arrKepETM = [arrKepETM; mapKep(maskCt)];
    arrVeETM = [arrVeETM; mapVe(maskCt)];
    arrVpETM = [arrVpETM; mapVp(maskCt)];
    % From Reference Region Model + RRIFT
    arrKtRRIFT = [arrKtRRIFT; mapKtR(maskCt)];
    arrKepRRIFT = [arrKepRRIFT; mapKepR(maskCt)];
    arrVeRRIFT = [arrVeRRIFT; mapVeR(maskCt)];
    arrVpRRIFT = [arrVpRRIFT; mapVpR(maskCt)];
    % From Reference Region Model + RRIFT with PopAvg AIF
    arrKtRRIFTPop = [arrKtRRIFTPop; mapKtR(maskCt).*estKtRRPop./estKtRR];
    arrVeRRIFTPop = [arrVeRRIFTPop; mapVeR(maskCt).*estVeRRPop./estVeRR];
    arrVpRRIFTPop = [arrVpRRIFTPop; mapVpR(maskCt).*estKtRRPop./estKtRR];
    % From Reference Region Model (without RRIFT)
    arrKtRRM = [arrKtRRM; mapKtR(maskCt)./estKtRR];
    arrVeRRM = [arrVeRRM; mapVeR(maskCt)./estVeRR];
    arrVpRRM = [arrVpRRM; mapVpR(maskCt)./estKtRR];
%     % ^ These maps were multipled by estKtRR and estVeRR in a04_doRRIFT, so
%     % we have to undo the scaling to get the original fits without RRIFT
    %% Collect the curves for input function & muscle
    crrList = [crrList Crr];
    cpList = [cpList Cp];
    tList = [tList t];
    %% Compute the Concordance Correlation Coefficient for each dataset
    % This is the CCC for the tumour fits
    cccKt(i) = CCC(mapKt(maskCt),mapKtR(maskCt));
    cccKep(i) = CCC(mapKep(maskCt),mapKepR(maskCt));
    cccVe(i) = CCC(mapVe(maskCt),mapVeR(maskCt));
    cccVp(i) = CCC(mapVp(maskCt),mapVpR(maskCt));
    %% Concatenate the estimates for muscle (from RRIFT or Tofts model)
    % From Extended Tofts model
    arrKtETM_m(i) = ETM.muscle(1);
    arrKepETM_m(i) = ETM.muscle(2);
    arrVeETM_m(i) = ETM.muscle(1)./ETM.muscle(2);
    arrVpETM_m(i) = ETM.muscle(3);
    % From RRIFT
    arrKtRRIFT_m(i) = estKtRR;
    arrKepRRIFT_m(i) = estKepRR;
    arrVeRRIFT_m(i) = estVeRR;
    
    % From RRIFT
    arrKtRRIFTPop_m(i) = estKtRRPop;
    arrVeRRIFTPop_m(i) = estVeRRPop;
    
    % Confidence Intervals
    arrKtRRci(i) = ciKtRR;
    arrKepRRci(i) = ciKepRR;
    arrVeRRci(i) = ciVeRR;
    arrKtRRTci(i) = ciKtRRT;
    arrKepRRTci(i) = ciKepRRT;
    arrVeRRTci(i) = ciVeRRT;
    
    % From Differential versione of RRIFT
    arrKtRRIFTD_m(i) = estKtRRdiff;
    arrVeRRIFTD_m(i) = estVeRRdiff;
    % kepRRD and kepRRIFT will be the same, since kepRR is estimate by
    % reference region model and not by RRIFT
    %% Collect additional statistics
    nVoxGood(i) = numGoodVox;
    nVox(i) = numVox;
    nVoxRR(i) = sum(maskCrr(:)>0); % Number of voxels in muscle contour
    snrList(i) = snr;
    cnrList(i) = cnr;
    sigmaList(i) = sigmaCt;
    t1CpList(i) = T1Cp;
    t1CrrList(i) = T1Crr;
    t1CtList(i) = T1Ct;
end
%% Cleanup - Part 1: Re-organize some results into a structure for clarity

[tumour, muscle] = deal(struct());

% ETM parameters for tumour
tumour.Kt.ETM = arrKtETM;
tumour.Kep.ETM = arrKepETM;
tumour.Ve.ETM = arrVeETM;
tumour.Vp.ETM = arrVpETM;

% RRIFT parameters for tumour
tumour.Kt.RRIFT = arrKtRRIFT;
tumour.Kep.RRIFT = arrKepRRIFT;
tumour.Ve.RRIFT = arrVeRRIFT;
tumour.Vp.RRIFT = arrVpRRIFT;

% RRIFT (assumed tail) parameters for tumour
tumour.Kt.RRIFT_Pop = arrKtRRIFTPop;
tumour.Kep.RRIFT_Pop = arrKepRRIFT; % kep is unaffected by tail choice
tumour.Ve.RRIFT_Pop = arrVeRRIFTPop;
tumour.Vp.RRIFT_Pop = arrVpRRIFTPop;

% Relative parameters from RRM
tumour.Kt.RRM = arrKtRRM;
tumour.Kep.RRM = arrKepRRIFT; % kep is unaffected by RRIFT
tumour.Ve.RRM = arrVeRRM;
tumour.Vp.RRM = arrVpRRM;

% ETM parameters for muscle
muscle.Kt.ETM = arrKtETM_m;
muscle.Kep.ETM = arrKepETM_m;
muscle.Ve.ETM = arrVeETM_m;

% RRIFT parameters for muscle
muscle.Kt.RRIFT = arrKtRRIFT_m;
muscle.Kep.RRIFT = arrKepRRIFT_m;
muscle.Ve.RRIFT = arrVeRRIFT_m;

% RRIFT (assumed tail) parameters for muscle
muscle.Kt.RRIFT_Pop = arrKtRRIFTPop_m;
muscle.Kep.RRIFT_Pop = arrKepRRIFT_m;
muscle.Ve.RRIFT_Pop = arrVeRRIFTPop_m;

%% Cleanup - Part 2: Clear all unnecessary variables
clearvars -except tumour muscle ccc* nVox* *List *ci

%% FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig. 8 Scatter of muscle parameters - RRIFT vs Tofts Model
cccs = round([ ...
        CCC(muscle.Kep.ETM,muscle.Kep.RRIFT),...
        CCC(muscle.Kt.ETM,muscle.Kt.RRIFT),...
        CCC(muscle.Ve.ETM,muscle.Ve.RRIFT)...
        ], 3);

figure('Position',[300,300,1500,400])

sz = 80;

subplot(1,3,1)
hold on;
% x and y values
xv = muscle.Kep.ETM;
yv = muscle.Kep.RRIFT;
% error bar from confidence interval
xr = arrKepRRTci;
yr = arrKepRRci;
xr = repmat(xr/2,[1 2]);
yr = repmat(yr/2,[1 2]);

plot([0 1],[0 1],'--','LineWidth',2,'Color',[0.4,0.4,0.4]);
errorbar(xv,yv,yr(:,1),yr(:,2),xr(:,1),xr(:,2),'k','LineStyle','none','LineWidth', 2)
scatter(xv,yv,sz,'k','LineWidth',2,'MarkerFaceColor',[1 1 1])
pbaspect([1 1 1])
ylim([0.1 0.8])
xlim([0.1 0.8])
customizeFig(gca);
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh, ['CCC: ' num2str(cccs(1))],'location','southeast')
legend boxoff
title('kep_{RR}')
xlabel('Tofts')
ylabel('RRIFT')

subplot(1,3,2)
hold on;
xv = muscle.Kt.ETM;
yv = muscle.Kt.RRIFT;
xr = arrKtRRTci;
yr = arrKtRRci;
xr = repmat(xr/2,[1 2]);
yr = repmat(yr/2,[1 2]);

plot([0 0.3],[0 0.3],'--','LineWidth',2,'Color',[0.4,0.4,0.4]);
errorbar(xv,yv,yr(:,1),yr(:,2),xr(:,1),xr(:,2),'k','LineStyle','none','LineWidth', 2)
scatter(xv,yv,sz,'k','LineWidth',2,'MarkerFaceColor',[1 1 1])
pbaspect([1 1 1])
set(gca,'YTick',[0:0.1:0.3])
set(gca,'XTick',[0:0.1:0.3])
ylim([0 0.2])
xlim([0 0.2])
customizeFig(gca);
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh, ['CCC: ' num2str(cccs(2))],'location','southeast')
legend boxoff
title('Ktrans_{RR}')
xlabel('Tofts')
ylabel('RRIFT')

subplot(1,3,3)
hold on;
xv = muscle.Ve.ETM;
yv = muscle.Ve.RRIFT;
xr = arrVeRRTci;
yr = arrVeRRci;
xr = repmat(xr/2,[1 2]);
yr = repmat(yr/2,[1 2]);

plot([0 1],[0 1],'--','LineWidth',2,'Color',[0.4,0.4,0.4]);
errorbar(xv,yv,yr(:,1),yr(:,2),xr(:,1),xr(:,2),'k','LineStyle','none','LineWidth', 2)
scatter(xv,yv,sz,'k','LineWidth',2,'MarkerFaceColor',[1 1 1])
pbaspect([1 1 1])
set(gca,'YTick',[0:0.1:0.3])
set(gca,'XTick',[0:0.1:0.3])
ylim([0.1 0.3])
xlim([0.1 0.3])
customizeFig(gca);
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh, ['CCC: ' num2str(cccs(3))],'location','southeast')
legend boxoff
title('ve_{RR}')
xlabel('Tofts')
ylabel('RRIFT')

%% Mean of Muscle parameters
disp('------ Muscle Parameters ------')
disp('Mean and StdDev of KtransRR from Tofts model:')
disp([mean(muscle.Kt.ETM) std(muscle.Kt.ETM)])

disp('Mean and StdDev of KtransRR from RRIFT:')
disp([mean(muscle.Kt.RRIFT) std(muscle.Kt.RRIFT)])

disp('Mean and StdDev of KtransRR from RRIFT with PopAvg AIF:')
disp([mean(muscle.Kt.RRIFT_Pop) std(muscle.Kt.RRIFT_Pop)])

disp('Mean and StdDev of veRR from Tofts Model:')
disp([mean(muscle.Ve.ETM) std(muscle.Ve.ETM)])

disp('Mean and StdDev of veRR from RRIFT:')
disp([mean(muscle.Ve.RRIFT) std(muscle.Ve.RRIFT)])

disp('Mean and StdDev of veRR from RRIFT with PopAvg AIF:')
disp([mean(muscle.Ve.RRIFT_Pop) std(muscle.Ve.RRIFT_Pop)])

%% Concordance Correlation Coefficient (CCC) for Reference Tissue
disp('------ Cocordance Correlation Coefficient for Muscle ------')
% Comparing estimates for reference tussye from Tofts Model Fit vs RRIFT
disp('CCC for KepRR, KtransRR, and veRR - between Tofts and RRIFT')
disp([CCC(muscle.Kep.ETM,muscle.Kep.RRIFT) CCC(muscle.Kt.ETM,muscle.Kt.RRIFT) CCC(muscle.Ve.ETM,muscle.Ve.RRIFT)])
disp('CCC for KepRR, KtransRR, and veRR - between Tofts and RRIFT with PopAvg AIF')
disp([CCC(muscle.Kep.ETM,muscle.Kep.RRIFT_Pop) CCC(muscle.Kt.ETM,muscle.Kt.RRIFT_Pop) CCC(muscle.Ve.ETM,muscle.Ve.RRIFT_Pop)])

%% Plot curves of input function and muscle for all patients
figure('Position',[300,300,1500,450])
subplot(1,2,1)
plot(tList,cpList,'LineWidth',1.5);
hold on;
popT = 1:1:max(tList(:))*60;
popT = popT'/60;
[Cp, Cb] = GeorgiouAif(popT,tList(7,1));
plot(popT,Cp,'--k','LineWidth',3);
ylim([-0.05 8])
ylabel('Concentration [mM]')
xlabel('Time [min]')
customizeFig(gca);

subplot(1,2,2)
plot(tList,crrList,'LineWidth',1.5)
ylabel('Concentration [mM]')
xlabel('Time [min]')
customizeFig(gca);

%% [Figs. 5 & S2] 2D Histograms - RRIFT vs Tofts Model
overlayInfo = [];

doLog = 0;


figure('Position',[100,300,1800,400]);

% 2D hist kt

valsA = tumour.Kt.ETM;
valsB = tumour.Kt.RRIFT;

valsA(valsA<0) = NaN;
valsB(valsB<0) = NaN;

if doLog
    valsA = log10(valsA);
    valsB = log10(valsB);
    valsA(valsA < -5) = NaN;
    valsB(valsB < -5) = NaN;
else
    valsA(valsA > 0.2) = NaN;
    valsB(valsB > 0.2) = NaN;
end

subplot(1,3,1)
h=histogram2(valsA,valsB,100,'DisplayStyle','tile','ShowEmptyBins','on');
[a1,a2] = CCC(valsA,valsB);
overlayInfo(1,1) = a1;
overlayInfo(2,1) = a2;
sum(isfinite(h.Data(:)))/2./length(h.Data);
imagesc(h.XBinEdges,h.YBinEdges,imgaussfilt(log10(h.Values'+1),0.5))
set(gca,'YDir','normal')
colormap('jet')
caxis([0 3])
hold on;
plot([min(valsA) max(valsA)],[min(valsA) max(valsA)],'w')
xlabel('ExtToftsModel')
ylabel('RRIFT (Measured Tail)')
title('kt')
colorbar
pbaspect([1 1 1])
%set(gca,'YTick',0:0.05:0.2)
customizeFig(gca);

% 2D hist ve
valsA = tumour.Ve.ETM;
valsB = tumour.Ve.RRIFT;

valsA(valsA<0) = NaN;
valsB(valsB<0) = NaN;

if doLog
    valsA = log10(valsA);
    valsB = log10(valsB);
    valsA(valsA < -5) = NaN;
    valsB(valsB < -5) = NaN;
else
    valsA(valsA > 0.5) = NaN;
    valsB(valsB > 0.5) = NaN;
end

subplot(1,3,2)
h=histogram2(valsA,valsB,100,'DisplayStyle','tile','ShowEmptyBins','on');
imagesc(h.XBinEdges,h.YBinEdges,imgaussfilt(log10(h.Values'+1),0.5))
[a1,a2] = CCC(valsA,valsB);
overlayInfo(1,2) = a1;
overlayInfo(2,2) = a2;
set(gca,'YDir','normal')
colormap('jet')
caxis([0 3])
hold on;
plot([min(valsA) max(valsA)],[min(valsA) max(valsA)],'w')
xlabel('ExtToftsModel')
ylabel('RRIFT (Measured Tail)')
title('ve')
colorbar
pbaspect([1 1 1])
%set(gca,'YTick',0:0.1:0.5)
customizeFig(gca);

% 2D hist vp

valsA = tumour.Vp.ETM;
valsB = tumour.Vp.RRIFT;

valsA(valsA<0) = NaN;
valsB(valsB<0) = NaN;

if doLog
    valsA = log10(valsA);
    valsB = log10(valsB);
    valsA(valsA < -5) = NaN;
    valsB(valsB < -5) = NaN;
else
    valsA(valsA > 0.05) = NaN;
    valsB(valsB > 0.05) = NaN;
end

subplot(1,3,3)
h=histogram2(valsA,valsB,100,'DisplayStyle','tile','ShowEmptyBins','on');
[a1,a2] = CCC(valsA,valsB);
overlayInfo(1,3) = a1;
overlayInfo(2,3) = a2;
sum(isfinite(h.Data(:)))/2./length(h.Data);
imagesc(h.XBinEdges,h.YBinEdges,imgaussfilt(log10(h.Values'+1),0.5))
set(gca,'YDir','normal')
colormap('jet')
caxis([0 3])
hold on;
plot([min(valsA) max(valsA)],[min(valsA) max(valsA)],'w')
xlabel('ExtToftsModel')
ylabel('RRIFT (Measured Tail)')
title('vp')
colorbar
pbaspect([1 1 1])
%set(gca,'YTick',0:0.01:0.05)
customizeFig(gca);

disp('------------')
disp('RRIFT with measured AIF-tail.') 
disp('CCC and Pearson correlation coefficient for Ktrans, ve, and vp')
disp(overlayInfo)
disp('')
%% [Figs. 5 & S2] 2D Histogram - RRIFT w/ PopAvg AIF vs Tofts Model
overlayInfo = [];

doLog = 0;


figure('Position',[100,300,1800,400]);

% 2D hist kt

valsA = tumour.Kt.ETM;
valsB = tumour.Kt.RRIFT_Pop;

valsA(valsA<0) = NaN;
valsB(valsB<0) = NaN;

if doLog
    valsA = log10(valsA);
    valsB = log10(valsB);
    valsA(valsA < -5) = NaN;
    valsB(valsB < -5) = NaN;
else
    valsA(valsA > 0.2) = NaN;
    valsB(valsB > 0.2) = NaN;
end

subplot(1,3,1)
h=histogram2(valsA,valsB,100,'DisplayStyle','tile','ShowEmptyBins','on');
[a1,a2] = CCC(valsA,valsB);
overlayInfo(1,1) = a1;
overlayInfo(2,1) = a2;
sum(isfinite(h.Data(:)))/2./length(h.Data);
imagesc(h.XBinEdges,h.YBinEdges,imgaussfilt(log10(h.Values'+1),0.5))
set(gca,'YDir','normal')
colormap('jet')
caxis([0 3])
hold on;
plot([min(valsA) max(valsA)],[min(valsA) max(valsA)],'w')
xlabel('ExtToftsModel')
ylabel('RRIFT (Assumed Tail)')
title('kt')
colorbar
pbaspect([1 1 1])
%set(gca,'YTick',0:0.05:0.2)
customizeFig(gca);

% 2D hist ve
valsA = tumour.Ve.ETM;
valsB = tumour.Ve.RRIFT_Pop;

valsA(valsA<0) = NaN;
valsB(valsB<0) = NaN;

if doLog
    valsA = log10(valsA);
    valsB = log10(valsB);
    valsA(valsA < -5) = NaN;
    valsB(valsB < -5) = NaN;
else
    valsA(valsA > 0.5) = NaN;
    valsB(valsB > 0.5) = NaN;
end

subplot(1,3,2)
h=histogram2(valsA,valsB,100,'DisplayStyle','tile','ShowEmptyBins','on');
imagesc(h.XBinEdges,h.YBinEdges,imgaussfilt(log10(h.Values'+1),0.5))
[a1,a2] = CCC(valsA,valsB);
overlayInfo(1,2) = a1;
overlayInfo(2,2) = a2;
set(gca,'YDir','normal')
colormap('jet')
caxis([0 3])
hold on;
plot([min(valsA) max(valsA)],[min(valsA) max(valsA)],'w')
xlabel('ExtToftsModel')
ylabel('RRIFT (Assumed Tail)')
title('ve')
colorbar
pbaspect([1 1 1])
%set(gca,'YTick',0:0.1:0.5)
customizeFig(gca);

% 2D hist vp

valsA = tumour.Vp.ETM;
valsB = tumour.Vp.RRIFT_Pop;

valsA(valsA<0) = NaN;
valsB(valsB<0) = NaN;

if doLog
    valsA = log10(valsA);
    valsB = log10(valsB);
    valsA(valsA < -5) = NaN;
    valsB(valsB < -5) = NaN;
else
    valsA(valsA > 0.05) = NaN;
    valsB(valsB > 0.05) = NaN;
end

subplot(1,3,3)
h=histogram2(valsA,valsB,100,'DisplayStyle','tile','ShowEmptyBins','on');
[a1,a2] = CCC(valsA,valsB);
overlayInfo(1,3) = a1;
overlayInfo(2,3) = a2;
sum(isfinite(h.Data(:)))/2./length(h.Data);
imagesc(h.XBinEdges,h.YBinEdges,imgaussfilt(log10(h.Values'+1),0.5))
set(gca,'YDir','normal')
colormap('jet')
caxis([0 3])
hold on;
plot([min(valsA) max(valsA)],[min(valsA) max(valsA)],'w')
xlabel('ExtToftsModel')
ylabel('RRIFT (Assumed Tail)')
title('vp')
colorbar
pbaspect([1 1 1])
%set(gca,'YTick',0:0.01:0.05)
customizeFig(gca);

disp('RRIFT with assumed AIF-tail.') 
disp('CCC and Pearson correlation coefficient for Ktrans, ve, and vp')
disp(overlayInfo)
disp('')

%% Histogram - assumed values of 0.07/0.14 (same as simulation)
overlayInfo = [];
doLog = 0;

assumedKtRR = 0.07;
assumedVeRR = 0.14;

figure('Position',[100,300,1800,400]);

% 2D hist kt

valsA = tumour.Kt.ETM;
valsB = tumour.Kt.RRM*0.07;

valsA(valsA<0) = NaN;
valsB(valsB<0) = NaN;

if doLog
    valsA = log10(valsA);
    valsB = log10(valsB);
    valsA(valsA < -5) = NaN;
    valsB(valsB < -5) = NaN;
else
    valsA(valsA > 0.2) = NaN;
    valsB(valsB > 0.2) = NaN;
end

subplot(1,3,1)
h=histogram2(valsA,valsB,100,'DisplayStyle','tile','ShowEmptyBins','on');
[a1,a2] = CCC(valsA,valsB);
overlayInfo(1,1) = a1;
overlayInfo(2,1) = a2;
sum(isfinite(h.Data(:)))/2./length(h.Data);
imagesc(h.XBinEdges,h.YBinEdges,imgaussfilt(log10(h.Values'+1),0.5))
set(gca,'YDir','normal')
colormap('jet')
caxis([0 3])
hold on;
plot([min(valsA) max(valsA)],[min(valsB) max(valsB)],'w')
xlabel('ExtToftsModel')
ylabel('RRM (Assumed KtransRR=0.07)')
title('kt')
colorbar
pbaspect([1 1 1])
%set(gca,'YTick',0:0.05:0.2)

% 2D hist ve
valsA = tumour.Ve.ETM;
valsB = tumour.Ve.RRM*0.14;

valsA(valsA<0) = NaN;
valsB(valsB<0) = NaN;

if doLog
    valsA = log10(valsA);
    valsB = log10(valsB);
    valsA(valsA < -5) = NaN;
    valsB(valsB < -5) = NaN;
else
    valsA(valsA > 0.5) = NaN;
    valsB(valsB > 0.5) = NaN;
end

subplot(1,3,2)
h=histogram2(valsA,valsB,100,'DisplayStyle','tile','ShowEmptyBins','on');
imagesc(h.XBinEdges,h.YBinEdges,imgaussfilt(log10(h.Values'+1),0.5))
[a1,a2] = CCC(valsA,valsB);
overlayInfo(1,2) = a1;
overlayInfo(2,2) = a2;
set(gca,'YDir','normal')
colormap('jet')
caxis([0 3])
hold on;
plot([min(valsA) max(valsA)],[min(valsA) max(valsA)],'w')
xlabel('ExtToftsModel')
ylabel('RRM (Assumed veRR=0.14)')
title('ve')
colorbar
pbaspect([1 1 1])
%set(gca,'YTick',0:0.1:0.5)

% 2D hist vp

valsA = tumour.Vp.ETM;
valsB = tumour.Vp.RRM*0.07;

valsA(valsA<0) = NaN;
valsB(valsB<0) = NaN;

if doLog
    valsA = log10(valsA);
    valsB = log10(valsB);
    valsA(valsA < -5) = NaN;
    valsB(valsB < -5) = NaN;
else
    valsA(valsA > 0.05) = NaN;
    valsB(valsB > 0.05) = NaN;
end

subplot(1,3,3)
h=histogram2(valsA,valsB,100,'DisplayStyle','tile','ShowEmptyBins','on');
[a1,a2] = CCC(valsA,valsB);
overlayInfo(1,3) = a1;
overlayInfo(2,3) = a2;
sum(isfinite(h.Data(:)))/2./length(h.Data);
imagesc(h.XBinEdges,h.YBinEdges,imgaussfilt(log10(h.Values'+1),0.5))
set(gca,'YDir','normal')
colormap('jet')
caxis([0 3])
hold on;
plot([min(valsA) max(valsA)],[min(valsA) max(valsA)],'w')
xlabel('ExtToftsModel')
ylabel('RRM (Assumed KtransRR=0.07)')
title('vp')
colorbar
pbaspect([1 1 1])
%set(gca,'YTick',0:0.01:0.05)

disp('RRM with assumed muscle parameters.') 
disp('CCC and Pearson correlation coefficient for Ktrans, ve, and vp')
disp(overlayInfo)
disp('')
%%
disp('Sigma (Standard deviation of noise in mM)')
disp([mean(sigmaList) std(sigmaList)])

disp('Contrast Noise Ratio (CNR)')
disp([mean(cnrList) std(cnrList)])

disp('Signal Noise Ratio (SNR)')
disp([mean(snrList) std(snrList)])
