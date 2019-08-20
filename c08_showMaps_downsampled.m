% Shows the maps for Ktrans, ve, vp from RRIFT and ETM at a high and low
% temporal resolution

clearvars
fclose('all')
addpath('./mfiles')

chosenIdx = 7;

inDir = './data/TCGA-GBM-Results/c06_downsampled';


matFiles = dir([inDir '/*.mat']);

for q=1:length(matFiles)
    if q~=chosenIdx
        continue
    end
    
    curFile = matFiles(q).name;
    load(fullfile(inDir,curFile));
    %%
    pkCE1 = pkCE(:,:,1);
    pkCE2 = pkCE(:,:,6);
    pkET1 = pkETM(:,:,1);
    pkET2 = pkETM(:,:,6);
    
    estVe1 = squeeze(pkET1(:,1)./pkET1(:,2));
    estVe2 = squeeze(pkET2(:,1)./pkET2(:,2));
    
    [sX, sY, sZ] = size(maskCt);
    
    [mapKtE1, mapKtE2, mapVeE1, mapVeE2, mapVpE1, mapVpE2,...
        mapKtR1, mapKtR2, mapVeR1, mapVeR2, mapVpR1, mapVpR2] = deal(zeros(sX,sY,sZ));
    
    mapKtE1(maskCt) = pkET1(:,1);
    mapKtE2(maskCt) = pkET2(:,1);
    
    mapVeE1(maskCt) = estVe1;
    mapVeE2(maskCt) = estVe2;
    
    mapVpE1(maskCt) = pkET1(:,3);
    mapVpE2(maskCt) = pkET2(:,3);
    
    mapKtR1(maskCt) = pkCE1(:,1);
    mapKtR2(maskCt) = pkCE2(:,1);
    
    mapVeR1(maskCt) = pkCE1(:,2);
    mapVeR2(maskCt) = pkCE2(:,2);
    
    mapVpR1(maskCt) = pkCE1(:,4);
    mapVpR2(maskCt) = pkCE2(:,4);
    
    %%
    mapKtE1 = AutoCrop(mapKtE1);
    mapKtE2 = AutoCrop(mapKtE2);
    mapVeE1 = AutoCrop(mapVeE1);
    mapVeE2 = AutoCrop(mapVeE2);
    mapVpE1 = AutoCrop(mapVpE1);
    mapVpE2 = AutoCrop(mapVpE2);
    
    mapKtR1 = AutoCrop(mapKtR1);
    mapKtR2 = AutoCrop(mapKtR2);
    mapVeR1 = AutoCrop(mapVeR1);
    mapVeR2 = AutoCrop(mapVeR2);
    mapVpR1 = AutoCrop(mapVpR1);
    mapVpR2 = AutoCrop(mapVpR2);
    %%
    mapVeE1(mapVeE1>1) = NaN;
    mapVeE2(mapVeE2>1) = NaN;
    
    mapVeR1(mapVeR1>1) = NaN;
    mapVeR2(mapVeR2>1) = NaN;
    
    mapVeE1(mapVeE1<0) = NaN;
    mapVeE2(mapVeE2<0) = NaN;
    
    mapVeR1(mapVeR1<0) = NaN;
    mapVeR2(mapVeR2<0) = NaN;
    %%
    mapVeE1(~isfinite(mapVeE1))=0;
	mapVeE2(~isfinite(mapVeE2))=0;
    mapVeR1(~isfinite(mapVeR1))=0;
	mapVeR2(~isfinite(mapVeR2))=0;
    %% Apply small amount of smoothing to suppress noisy fluctuations
    h=fspecial('disk',1);
    
    mapKtE1 = imfilter(mapKtE1,h);
    mapKtE2 = imfilter(mapKtE2,h);
    mapVeE1 = imfilter(mapVeE1,h);
    mapVeE2 = imfilter(mapVeE2,h);
    mapVpE1 = imfilter(mapVpE1,h);
    mapVpE2 = imfilter(mapVpE2,h);
    
    mapKtR1 = imfilter(mapKtR1,h);
    mapKtR2 = imfilter(mapKtR2,h);
    mapVeR1 = imfilter(mapVeR1,h);
    mapVeR2 = imfilter(mapVeR2,h);
    mapVpR1 = imfilter(mapVpR1,h);
    mapVpR2 = imfilter(mapVpR2,h);
end
%% Show maps for one slice
myS = 4;
cMax = 0.15;

figure

subplot(2,2,1)
imagesc(mapKtE1(:,:,myS))
caxis([0 cMax]);
axis image
axis off
colormap('jet')
title('ETM-Ktrans-Hi')

subplot(2,2,2)
imagesc(mapKtR1(:,:,myS))
caxis([0 cMax]);
axis image
axis off
title('RRIFT-Ktrans-Hi')

subplot(2,2,3)
imagesc(mapKtE2(:,:,myS))
caxis([0 cMax]);
axis image
axis off
title('ETM-Ktrans-Lo')

subplot(2,2,4)
imagesc(mapKtR2(:,:,myS))
caxis([0 cMax]);
axis image
axis off
title('RRIFT-Ktrans-Lo')

figure

cMax = .5;
subplot(2,2,1)
imagesc(mapVeE1(:,:,myS))
caxis([0 cMax]);
axis image
axis off
colormap('jet')
title('ETM-ve-Hi')
subplot(2,2,2)
imagesc(mapVeR1(:,:,myS))
caxis([0 cMax]);
axis image
axis off
title('RRIFT-ve-Hi')
subplot(2,2,3)
imagesc(mapVeE2(:,:,myS))
caxis([0 cMax]);
axis image
axis off
title('ETM-ve-Lo')
subplot(2,2,4)
imagesc(mapVeR2(:,:,myS))
caxis([0 cMax]);
axis image
axis off
title('RRIFT-ve-Lo')

figure

cMax = 0.15;
subplot(2,2,1)
imagesc(mapVpE1(:,:,myS))
caxis([0 cMax]);
axis image
axis off
colormap('jet')
title('ETM-vp-Hi')
subplot(2,2,2)
imagesc(mapVpR1(:,:,myS))
caxis([0 cMax]);
axis image
axis off
title('RRIFT-vp-Hi')
subplot(2,2,3)
imagesc(mapVpE2(:,:,myS))
caxis([0 cMax]);
axis image
axis off
title('ETM-vp-Lo')
subplot(2,2,4)
imagesc(mapVpR2(:,:,myS))
caxis([0 cMax]);
axis image
axis off
title('RRIFT-vp-Lo')
