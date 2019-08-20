function [Ct, T1, M0, dceData, t, dceHdr, vfaHdr, flipAngles, flipData] = SortItOutForMe(dcePath, vfaPath)

%%
[dceData, dceHdr, dceAcqTime] = getDicomImage(dcePath); 
[flipData, vfaHdr, ~, flipAngles] = getDicomImage(vfaPath);
%%
[sX, sY, sZF] = size(flipData);
sF = length(unique(flipAngles));
sZ = sZF/sF;
if round(sZ)~=sZ
    estSZ = sum(flipAngles==5);
    flipData = flipData(:,:,estSZ+1:end);
    [sX, sY, sZF] = size(flipData);
    sZ = sZF/sF;
    flipAngles = flipAngles(estSZ+1:end);
    if round(sZ)~=sZ && unique(flipAngles)~=sF
        error('Mismatch between number of slices in flip angle data')
    end
end
%%
flipData = reshape(flipData,[sX sY sZ sF]);
flipAngles = reshape(flipAngles,[sZ sF]);
%%
[~,~,sZT] = size(dceData);
sT = sZT/sZ;
if round(sT)~=sT
    error('Mismatch between number of slices in DCE data')
end
dceData = reshape(dceData,[sX sY sZ sT]);
%%
timeStep = dceAcqTime(sZ+1) - dceAcqTime(1);
t=timeStep:timeStep:timeStep*sT;
t=t'/60/1000;
%%
tmpData = sum(dceData,4);
mask = tmpData > 0.05*max(tmpData(:));
mask = repmat(sum(mask,3),[1 1 sZ])>0;
clearvars tmpData
%%
[T1, M0, ~] = doDESPOT1(flipData,deg2rad(flipAngles(1,:)),vfaHdr.RepetitionTime,mask);
%%
if contains(dceHdr.ContrastBolusAgent,'magnevist')
    % Relaxivity of Magnevist (Gd-DTPA) from Pintaske et al. 2006
    r1=3.3;
else
    error('Could not interpret ContrastBolusAgent in DCE header');
end
numPrecontrastFrames=3;
Ct = Signal2Conc2D(unravel(dceData)',T1(:),dceHdr.RepetitionTime,dceHdr.FlipAngle,r1,numPrecontrastFrames);
Ct(:,~mask(:))=0;
Ct = reshape(Ct',[sX sY sZ sT]);
Ct = SubtractBaseline(Ct,numPrecontrastFrames);
%%
end
