%% LOAD
addpath('./mfiles')
clearvars

load('./data/simMap.mat');
load('./data/simResultsTRes15-varKtRR.mat');

estRRM = zeros(9,3,100,repF);
estRRIFT = estRRM;
estETM = permute(pkETM,[1 3 2 4]);
estETM(:,2,:,:) = estETM(:,1,:,:)./estETM(:,2,:,:);

errRRM = zeros(9,3,100*repF);
errRRIFT = errRRM;
errETM = errRRM;

estVeRR = estKtRR ./ estKepRRS;

for i=1:9
        estRRM(i,1,:,:) = pkCERRM(i,:,1,:) * 0.07;
        estRRM(i,2,:,:) = pkCERRM(i,:,2,:) * 0.14;
        estRRM(i,3,:,:) = pkCERRM(i,:,4,:) * 0.07;
        
        for k=1:repF
            estRRIFT(i,1,:,k) = pkCERRM(i,:,1,k) * estKtRR(i,k);
            estRRIFT(i,2,:,k) = pkCERRM(i,:,2,k) * estVeRR(i,k);
            estRRIFT(i,3,:,k) = pkCERRM(i,:,4,k) * estKtRR(i,k);
        end
        
        errRRM(i,1,:) = reshape(PercentError(squeeze(estRRM(i,1,:,:)),trueKt(:)),1,[]);
        errRRM(i,2,:) = reshape(PercentError(squeeze(estRRM(i,2,:,:)),trueVe(:)),1,[]);
        errRRM(i,3,:) = reshape(PercentError(squeeze(estRRM(i,3,:,:)),trueVp(:)),1,[]);
        
        errRRIFT(i,1,:) = reshape(PercentError(squeeze(estRRIFT(i,1,:,:)),trueKt(:)),1,[]);
        errRRIFT(i,2,:) = reshape(PercentError(squeeze(estRRIFT(i,2,:,:)),trueVe(:)),1,[]);
        errRRIFT(i,3,:) = reshape(PercentError(squeeze(estRRIFT(i,3,:,:)),trueVp(:)),1,[]);
        
        errETM(i,1,:) = reshape(PercentError(squeeze(estETM(i,1,:,:)),trueKt(:)),1,[]);
        errETM(i,2,:) = reshape(PercentError(squeeze(estETM(i,2,:,:)),trueVe(:)),1,[]);
        errETM(i,3,:) = reshape(PercentError(squeeze(estETM(i,3,:,:)),trueVp(:)),1,[]);
end
trueKtRR0=trueKtRR; % Save to tmp variable because next loading step will overwrite it
load('./data/simResultsTRes15-varVeRR.mat');
estVeRR = estKtRR ./ estKepRRS;

for i=1:9
        estRRM(i,2,:,:) = pkCERRM(i,:,2,:) * 0.14;
        
        for k=1:repF
            estRRIFT(i,2,:,k) = pkCERRM(i,:,2,k) * estVeRR(i,k);
        end
        
        errRRM(i,2,:) = reshape(PercentError(squeeze(estRRM(i,2,:,:)),trueVe(:)),1,[]);
        errRRIFT(i,2,:) = reshape(PercentError(squeeze(estRRIFT(i,2,:,:)),trueVe(:)),1,[]);
        errETM(i,2,:) = reshape(PercentError(squeeze(estETM(i,2,:,:)),trueVe(:)),1,[]);
end
trueKtRR = trueKtRR0;
%% Error plots
% Err in Ktrans
cSize = 10;
yrange = [-100 100];
myColours = [0.895, 0.350, 0.280;...
             0.160, 0.140, 0.280;...
             0.750, 0.750, 0.750];

err1 = squeeze(errRRM(:,1,:));
err2 = squeeze(errRRIFT(:,1,:));
err3 = squeeze(errETM(:,1,:));

avg1 = nanmedian(err1,2);
avg2 = nanmedian(err2,2);
avg3 = nanmedian(err3,2);

ci1 = quantile(err1,[.25 .75],2);
ci2 = quantile(err2,[.25 .75],2);
ci3 = quantile(err3,[.25 .75],2);

figure('Position',[300,300,1600,400]); 
subplot(1,3,1)
hold on;
errorbar(trueKtRR(:,1),avg3,avg3-ci3(:,1),ci3(:,2)-avg3,'linewidth',2,'CapSize',cSize,'color',myColours(3,:));
errorbar(trueKtRR(:,1),avg2,avg2-ci2(:,1),ci2(:,2)-avg2,'linewidth',2,'CapSize',cSize,'color',myColours(2,:));
errorbar(trueKtRR(:,1),avg1,avg1-ci1(:,1),ci1(:,2)-avg1,'linewidth',2,'CapSize',cSize,'color',myColours(1,:));
hold off;
ylim(yrange);
xticks(trueKtRR(1:2:end,1));
ylabel('Percent Error in Ktrans')
xlabel('True KtransRR')

% Err in Ve

err1 = squeeze(errRRM(:,2,:));
err2 = squeeze(errRRIFT(:,2,:));
err3 = squeeze(errETM(:,2,:));

avg1 = nanmedian(err1,2);
avg2 = nanmedian(err2,2);
avg3 = nanmedian(err3,2);

ci1 = quantile(err1,[.25 .75],2);
ci2 = quantile(err2,[.25 .75],2);
ci3 = quantile(err3,[.25 .75],2);

cSize = 10;
subplot(1,3,2)
hold on;
errorbar(trueVeRR,avg3,avg3-ci3(:,1),ci3(:,2)-avg3,'linewidth',2,'CapSize',cSize,'color',myColours(3,:));
errorbar(trueVeRR,avg2,avg2-ci2(:,1),ci2(:,2)-avg2,'linewidth',2,'CapSize',cSize,'color',myColours(2,:));
errorbar(trueVeRR,avg1,avg1-ci1(:,1),ci1(:,2)-avg1,'linewidth',2,'CapSize',cSize,'color',myColours(1,:));
hold off;
ylim(yrange);
ylabel('Percent Error in ve')
xlabel('True veRR')
legend('ETM','RRIFT','RRM w/ fixed RR params')
legend boxoff

% Err in Vp
err1 = squeeze(errRRM(:,3,:));
err2 = squeeze(errRRIFT(:,3,:));
err3 = squeeze(errETM(:,3,:));

avg1 = nanmedian(err1,2);
avg2 = nanmedian(err2,2);
avg3 = nanmedian(err3,2);

ci1 = quantile(err1,[.25 .75],2);
ci2 = quantile(err2,[.25 .75],2);
ci3 = quantile(err3,[.25 .75],2);

cSize = 10;
subplot(1,3,3)
hold on;
errorbar(trueKtRR(:,1),avg3,avg3-ci3(:,1),ci3(:,2)-avg3,'linewidth',2,'CapSize',cSize,'color',myColours(3,:));
errorbar(trueKtRR(:,1),avg2,avg2-ci2(:,1),ci2(:,2)-avg2,'linewidth',2,'CapSize',cSize,'color',myColours(2,:));
errorbar(trueKtRR(:,1),avg1,avg1-ci1(:,1),ci1(:,2)-avg1,'linewidth',2,'CapSize',cSize,'color',myColours(1,:));
hold off;
ylim(yrange);
xticks(trueKtRR(1:2:end,1));
ylabel('Percent Error in vp')
xlabel('True KtransRR')