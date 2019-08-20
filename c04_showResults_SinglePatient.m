%% Description
% Loads maps from RRIFT and Extended Tofts Model
% Displays maps for a single patient at 3 slices
% Also plots the AIF/Muscle curve and RRIFT fit

% Estimated runtime: 1 second

%%
clearvars
fclose('all')
addpath('./mfiles')
tic
inDir = './data/TCGA-GBM-Results/c02_postprocessed';
matFiles = dir([inDir '/*.mat']);

idx = 7; % Chosen Patient (choices are between 1:9)
curFile = matFiles(idx).name;

load(fullfile(inDir,curFile));

% Remove unphysical values from ve and vp  
mapVe(mapVe>1)=NaN;
mapVe(mapVe<0)=NaN;
mapVeR(mapVeR>1)=NaN;
mapVeR(mapVeR<0)=NaN;

mapVp(mapVp>1)=NaN;
mapVp(mapVp<0)=NaN;
mapVpR(mapVpR>1)=NaN;
mapVpR(mapVpR<0)=NaN;

%% Crop the image

mapKt = AutoCrop(mapKt);
mapKep = AutoCrop(mapKep);
mapVe = AutoCrop(mapVe);
mapVp = AutoCrop(mapVp);

mapKtR = AutoCrop(mapKtR);
mapKepR = AutoCrop(mapKepR);
mapVeR = AutoCrop(mapVeR);
mapVpR = AutoCrop(mapVpR);

%% Plot slices for P7
slices = [7,9,11];
nS = length(slices);
majula = jet(100);
% Ktrans
figure
clims = [0 0.16];
for i = 1:length(slices)
    subplot(2,nS,i)
    imagesc(mapKt(:,:,slices(i)));
    caxis(clims); colormap(majula); set(gca,'XColor','none','YColor','none');
    title('Ktrans-ETM')

    subplot(2,nS,nS+i)
    imagesc(mapKtR(:,:,slices(i)));
    caxis(clims); colormap(majula); set(gca,'XColor','none','YColor','none');
    title('Ktrans-RRIFT')
end

% ve
h=fspecial('disk',1);
mapVeF = mapVe;
mapVeFR = mapVeR;
mapVeF(~isfinite(mapVeF))=0;
mapVeFR(~isfinite(mapVeFR))=0;
mapVeF = imfilter(mapVeF,h);
mapVeFR = imfilter(mapVeFR,h);

figure
clims = [0 0.6];
for i = 1:length(slices)
    subplot(2,nS,i)
    imagesc(mapVeF(:,:,slices(i)));
    caxis(clims); colormap(majula); set(gca,'XColor','none','YColor','none');
    title('ve-ETM')

    subplot(2,nS,nS+i)
    imagesc(mapVeFR(:,:,slices(i)));
    caxis(clims); colormap(majula); set(gca,'XColor','none','YColor','none');
    title('ve-RRIFT')
end

% vp
figure
clims = [0 0.1];
for i = 1:length(slices)
    subplot(2,nS,i)
    imagesc(mapVp(:,:,slices(i)));
    caxis(clims); colormap(majula); set(gca,'XColor','none','YColor','none');
    title('vp-ETM')

    subplot(2,nS,nS+i)
    imagesc(mapVpR(:,:,slices(i)));
    caxis(clims); colormap(majula); set(gca,'XColor','none','YColor','none');
    title('vp-RRIFT')
end
%% Plot the Cp, Crr, and RRIFT fit for the patient
figure('Position',[300,300,1500,400])

subplot(1,3,1)
plot(t,Cp,'r','LineWidth',2)
customizeFig(gca);
xlim([0 7])
xlabel('Time [min]');
ylabel('Concentration [mM]')
title('Concentration in Blood Plasma')

subplot(1,3,2)
plot(t,Crr,'g','LineWidth',2)
xlim([0 7])
ylim([0 0.3])
customizeFig(gca);
xlabel('Time [min]');
ylabel('Concentration [mM]')
title('Concentration in Reference Region (Muscle)')

subplot(1,3,3)
scatter(denum,num,70,'r','LineWidth',2,'MarkerEdgeAlpha',.6)
hold on;
plot(denum,(denum\num)*denum,'-.k','LineWidth',2)
customizeFig(gca);
xlabel('Denominator')
ylabel('Numerator')
title('Example of RRIFT fit')
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
dummyh2 = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh, ['R^2: ' num2str(corr2(num, denum)^2)],'location','southeast')
legend boxoff
%% Display the estimates muscle parameters from RRIFT and from ETM fits
disp('Estimated muscle parameters from RRIFT [KtransRR veRR]:')
disp([estKtRR estVeRR])
disp('Estimates muscle parameters from ETM [KtransRR veRR vpRR]:')
disp([ETM.muscle(1) ETM.muscle(1)./ETM.muscle(2) ETM.muscle(3)]);
toc