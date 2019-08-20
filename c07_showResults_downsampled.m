% Plots the percent changes in Ktrans, ve, vp for RRIFT and ETM on
% downsampled in-vivo data.

% Estimated runtime: 2 seconds

clearvars
fclose('all')
addpath('./mfiles')

inDir = './data/TCGA-GBM-Results/c06_downsampled';


matFiles = dir([inDir '/*.mat']);

pkCEs = [];
pkETs = [];
estKepRR = [];
estKtRR = [];

for q=1:length(matFiles)
    curFile = matFiles(q).name;
    load(fullfile(inDir,curFile));
    %%
    pkCEs = [pkCEs; pkCE];
    pkETs = [pkETs; pkETM];
    estKepRR = [estKepRR; estKepRRs];
    estKtRR = [estKtRR; estKtRRs];
end
%%
errKtCE = PercentError(squeeze(pkCEs(:,1,:)),pkCEs(:,1,1));
errVeCE = PercentError(squeeze(pkCEs(:,2,:)),pkCEs(:,2,1));
errKepCE = PercentError(squeeze(pkCEs(:,3,:)),pkCEs(:,3,1));
errVpCE = PercentError(squeeze(pkCEs(:,4,:)),pkCEs(:,4,1));

estVe = squeeze(pkETs(:,1,:)./pkETs(:,2,:));

errKtET = PercentError(squeeze(pkETs(:,1,:)),pkETs(:,1,1));
errVeET = PercentError(estVe(:,:),estVe(:,1));
errKepET = PercentError(squeeze(pkETs(:,2,:)),pkETs(:,2,1));
errVpET = PercentError(squeeze(pkETs(:,3,:)),pkETs(:,3,1));
%%
cSize = 10;

% Err Kt
figure('Position',[300,300,1500,400]); 

curErr = errKtCE;
errQt = quantile(curErr,[.25 .75]);
errMd = nanmedian(curErr);

subplot(1,3,1)
errorbar(TRes,errMd,abs(errQt(1,:)-errMd),abs(errQt(2,:)-errMd),'linewidth',2,'CapSize',cSize);

curErr = errKtET;
errQt = quantile(curErr,[.25 .75]);
errMd = nanmedian(curErr);

hold on;
errorbar(TRes,errMd,abs(errQt(1,:)-errMd),abs(errQt(2,:)-errMd),'linewidth',2,'CapSize',cSize);

ylim([-100 100]);
title('Ktrans')
xlabel('Temporal Resolutions [s]')
ylabel('Percent Error')

% Err Ve
curErr = errVeCE;
errQt = quantile(curErr,[.25 .75]);
errMd = nanmedian(curErr);

subplot(1,3,2)
errorbar(TRes,errMd,abs(errQt(1,:)-errMd),abs(errQt(2,:)-errMd),'linewidth',2,'CapSize',cSize);

curErr = errVeET;
errQt = quantile(curErr,[.25 .75]);
errMd = nanmedian(curErr);

hold on;
errorbar(TRes,errMd,abs(errQt(1,:)-errMd),abs(errQt(2,:)-errMd),'linewidth',2,'CapSize',cSize);
ylim([-100 100]);

title('ve')
xlabel('Temporal Resolutions [s]')
legend('RRIFT','ETM')
legend boxoff

% Err Vp
curErr = errVpCE;
errQt = quantile(curErr,[.25 .75]);
errMd = nanmedian(curErr);

subplot(1,3,3)
errorbar(TRes,errMd,abs(errQt(1,:)-errMd),abs(errQt(2,:)-errMd),'linewidth',2,'CapSize',cSize);

curErr = errVpET;
errQt = quantile(curErr,[.25 .75]);
errMd = nanmedian(curErr);

hold on;
errorbar(TRes,errMd,abs(errQt(1,:)-errMd),abs(errQt(2,:)-errMd),'linewidth',2,'CapSize',cSize);
ylim([-100 100]);

title('vp')
xlabel('Temporal Resolutions [s]')