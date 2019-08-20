addpath('./mfiles')
clearvars

load('./data/simMap.mat');
load('./data/simResults.mat');
%%
[nVox nP nRep nRes nSig] = size(params.CERRM);

pkCE = reshape(params.CERRM,[nVox 5 nRep*nRes*nSig]);
pkET = reshape(params.ETM,[nVox 3 nRep*nRes*nSig]);
pkCE = permute(pkCE,[1 3 2]);
pkET = permute(pkET,[1 3 2]);
%%
rrKt = repmat(estKtRR(:),[1 nVox])';
rrKep = repmat(estKepRRS(:),[1 nVox])';
rrVe = rrKt./rrKep;
%%
errKtCE = PercentError(pkCE(:,:,1).*rrKt,trueKt(:));
errVeCE = PercentError(pkCE(:,:,2).*rrVe,trueVe(:));
errVpCE = PercentError(pkCE(:,:,4).*rrKt,trueVp(:));

errKtET = PercentError(pkET(:,:,1),trueKt(:));
errVeET = PercentError(pkET(:,:,1)./pkET(:,:,2),trueVe(:));
errVpET = PercentError(pkET(:,:,3),trueVp(:));
%%
errKtCE = reshape(errKtCE,[nVox*nRep nRes nSig]);
errVeCE = reshape(errVeCE,[nVox*nRep nRes nSig]);
errVpCE = reshape(errVpCE,[nVox*nRep nRes nSig]);

errKtET = reshape(errKtET,[nVox*nRep nRes nSig]);
errVeET = reshape(errVeET,[nVox*nRep nRes nSig]);
errVpET = reshape(errVpET,[nVox*nRep nRes nSig]);
%%
errKepRR = PercentError(reshape(estKepRRS,[nRep nRes*nSig]),kepRR);
errKtRR = PercentError(reshape(estKtRR,[nRep nRes*nSig]),ktRR);
errVeRR = PercentError(reshape(estKtRR./estKepRRS,[nRep nRes*nSig]),veRR);

errKepRRD = PercentError(reshape(estKepRRS,[nRep nRes*nSig]),kepRR);
errKtRRD = PercentError(reshape(estKtRRD,[nRep nRes*nSig]),ktRR);
errVeRRD = PercentError(reshape(estKtRRD./estKepRRS,[nRep nRes*nSig]),veRR);

errKepRR = reshape(errKepRR,[nRep nRes nSig]);
errKtRR = reshape(errKtRR,[nRep nRes nSig]);
errVeRR = reshape(errVeRR,[nRep nRes nSig]);

errKepRRD = reshape(errKepRRD,[nRep nRes nSig]);
errKtRRD = reshape(errKtRRD,[nRep nRes nSig]);
errVeRRD = reshape(errVeRRD,[nRep nRes nSig]);
%%
errKtCE = shiftdim(errKtCE,1);
errVeCE = shiftdim(errVeCE,1);
errVpCE = shiftdim(errVpCE,1);

errKtET = shiftdim(errKtET,1);
errVeET = shiftdim(errVeET,1);
errVpET = shiftdim(errVpET,1);

errKepRR = shiftdim(errKepRR,1);
errKtRR = shiftdim(errKtRR,1);
errVeRR = shiftdim(errVeRR,1);

errKepRRD = shiftdim(errKepRRD,1);
errKtRRD = shiftdim(errKtRRD,1);
errVeRRD = shiftdim(errVeRRD,1);
%% Median/IQR RR
cSize = 10;
iList = 1:4;

myColours = [0.000, 0.408, 0.216;...
             0.192, 0.640, 0.329;...
             0.471, 0.777, 0.474;...
             0.761, 0.902, 0.600]; 

figure('Position',[300,300,1600,400]); 

subplot(1,3,1); hold on;
curErr = errKepRR;
errQt = quantile(curErr,[.25 .75],3);
errMd = median(curErr,3);
for i=iList
errorbar(listSigmaC,errMd(i,:),abs(errQt(i,:,1)-errMd(i,:)),abs(errQt(i,:,2)-errMd(i,:)),'linewidth',2,'CapSize',cSize, 'color',myColours(i,:));
end
ylim([-35 5])
set(gca,'XTick',[0:0.01:0.05])
title('kepRR')
xlabel('Noise [mM]')
ylabel('Percent Error')
customizeFig(gca);

subplot(1,3,2); hold on;
curErr = errKtRR;
errQt = quantile(curErr,[.25 .75],3);
errMd = median(curErr,3);
for i=iList
errorbar(listSigmaC,errMd(i,:),abs(errQt(i,:,1)-errMd(i,:)),abs(errQt(i,:,2)-errMd(i,:)),'linewidth',2,'CapSize',cSize,'color',myColours(i,:));
end
ylim([-35 5])
set(gca,'XTick',[0:0.01:0.05])
title('ktransRR')
xlabel('Noise [mM]')
customizeFig(gca);

subplot(1,3,3); hold on;
curErr = errVeRR;
errQt = quantile(curErr,[.25 .75],3);
errMd = median(curErr,3);
for i=iList
errorbar(listSigmaC,errMd(i,:),abs(errQt(i,:,1)-errMd(i,:)),abs(errQt(i,:,2)-errMd(i,:)),'linewidth',2,'CapSize',cSize,'color',myColours(i,:));
end
ylim([-35 5])
set(gca,'XTick',[0:0.01:0.05])
title('veRR')
xlabel('Noise [mM]')
legend({num2str(TRes(iList(1))),num2str(TRes(iList(2))),num2str(TRes(iList(3))),num2str(TRes(iList(4)))},'location','southeast')
legend boxoff
customizeFig(gca);
%% Median/IQR for TOI
cSize = 10;
iList = 1:4;

myColours = [0.145, 0.204, 0.580;...
             0.172, 0.498, 0.722;...
             0.255, 0.714, 0.769;...
             0.631, 0.855, 0.706];

figure('Position',[300,300,1500,400]); 

subplot(1,3,1); hold on;
curErr = errKtCE;
errQt = quantile(curErr,[.25 .75],3);
errMd = median(curErr,3);
for i=iList
errorbar(listSigmaC,errMd(i,:),abs(errQt(i,:,1)-errMd(i,:)),abs(errQt(i,:,2)-errMd(i,:)),'linewidth',2,'CapSize',cSize,'color',myColours(i,:));
end
ylim([-50 50])
set(gca,'XTick',[0:0.01:0.05])
title('ktrans')
xlabel('Noise [mM]')
ylabel('Percent Error')
legend({num2str(TRes(iList(1))),num2str(TRes(iList(2))),num2str(TRes(iList(3))),num2str(TRes(iList(4)))},'location','southeast')
legend boxoff

subplot(1,3,2); hold on;
curErr = errVeCE;
errQt = quantile(curErr,[.25 .75],3);
errMd = median(curErr,3);
for i=iList
errorbar(listSigmaC,errMd(i,:),abs(errQt(i,:,1)-errMd(i,:)),abs(errQt(i,:,2)-errMd(i,:)),'linewidth',2,'CapSize',cSize,'color',myColours(i,:));
end
ylim([-50 50])
set(gca,'XTick',[0:0.01:0.05])
title('ve')
xlabel('Noise [mM]')


subplot(1,3,3); hold on;
curErr = errVpCE;
errQt = quantile(curErr,[.25 .75],3);
errMd = median(curErr,3);
for i=iList
errorbar(listSigmaC,errMd(i,:),abs(errQt(i,:,1)-errMd(i,:)),abs(errQt(i,:,2)-errMd(i,:)),'linewidth',2,'CapSize',cSize,'color',myColours(i,:));
end
ylim([-50 50])
set(gca,'XTick',[0:0.01:0.05])
title('vp')
xlabel('Noise [mM]')

%% IQR stats - RR

disp('Median err kepRR')
disp(median(errKepRR,3))
disp('Median err ktRR')
disp(median(errKtRR,3))
disp('Median err veRR')
disp(median(errVeRR,3))

disp('IQR err kepRR')
disp(iqr(errKepRR,3))
disp('IQR err ktRR')
disp(iqr(errKtRR,3))
disp('IQR err veRR')
disp(iqr(errVeRR,3))

%% IQR stats - TOI

disp('RRIFT ====================')
disp('Median err kt')
disp(median(errKtCE,3))
disp('Median err ve')
disp(median(errVeCE,3))
disp('Median err vp')
disp(median(errVpCE,3))

disp('IQR err kt')
disp(iqr(errKtCE,3))
disp('IQR err ve')
disp(iqr(errVeCE,3))
disp('IQR err vp')
disp(iqr(errVpCE,3))

%% IQR stats - TOI - ETM

disp('ETM ====================')
disp('Median err kt')
disp(median(errKtET,3))
disp('Median err ve')
disp(median(errVeET,3))
disp('Median err vp')
disp(median(errVpET,3))

disp('IQR err kt')
disp(iqr(errKtET,3))
disp('IQR err ve')
disp(iqr(errVeET,3))
disp('IQR err vp')
disp(iqr(errVpET,3))

%%
return
%% EXTRA - Run blocks with Ctrl+Enter
%% 

disp('RRIFT - ETM ====================')
disp('Median err kt')
disp(abs(median(errKtCE,3)) - abs(median(errKtET,3)))
disp('Median err ve')
disp(abs(median(errVeCE,3)) - abs(median(errVeET,3)))
disp('Median err vp')
disp(abs(median(errVpCE,3)) - abs(median(errVpET,3)))
%%
disp('IQR err kt')
disp(iqr(errKtET,3))
disp('IQR err ve')
disp(iqr(errVeET,3))
disp('IQR err vp')
disp(iqr(errVpET,3))
%% [EXTRA] Median/IQR for TOI using ETM
cSize = 10;
figure('Position',[300,300,1500,400]); 

subplot(1,3,1); hold on;
curErr = errKtET;
errQt = quantile(curErr,[.25 .75],3);
errMd = median(curErr,3);
for i=iList
errorbar(listSigmaC,errMd(i,:),abs(errQt(i,:,1)-errMd(i,:)),abs(errQt(i,:,2)-errMd(i,:)),'linewidth',2,'CapSize',cSize);
end
ylim([-50 50])
set(gca,'XTick',[0:0.01:0.05])
title('ktrans')
xlabel('Noise [mM]')
ylabel('Percent Error')

subplot(1,3,2); hold on;
curErr = errVeET;
errQt = quantile(curErr,[.25 .75],3);
errMd = median(curErr,3);
for i=iList
errorbar(listSigmaC,errMd(i,:),abs(errQt(i,:,1)-errMd(i,:)),abs(errQt(i,:,2)-errMd(i,:)),'linewidth',2,'CapSize',cSize);
end
ylim([-50 50])
set(gca,'XTick',[0:0.01:0.05])
title('ve')
xlabel('Noise [mM]')

subplot(1,3,3); hold on;
curErr = errVpET;
errQt = quantile(curErr,[.25 .75],3);
errMd = median(curErr,3);
for i=iList
errorbar(listSigmaC,errMd(i,:),abs(errQt(i,:,1)-errMd(i,:)),abs(errQt(i,:,2)-errMd(i,:)),'linewidth',2,'CapSize',cSize);
end
ylim([-50 50])
set(gca,'XTick',[0:0.01:0.05])
title('vp')
xlabel('Noise [mM]')
legend({num2str(TRes(iList(1))),num2str(TRes(iList(2))),num2str(TRes(iList(3))),num2str(TRes(iList(4)))},'location','southeast')
legend boxoff