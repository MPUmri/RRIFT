% This script does RRIFT on in-vivo data but changes the start of the tail.
% The goal is to see how changing the tail duration affects the KtransRR
% estimate.

% Estimated runtime: <30 seconds

%%
clearvars
fclose('all')
addpath('./mfiles')

inDir = './data/TCGA-GBM-Results/c01_preprocessed';

% There are 70 frames, so tail can begin from frames 1:69
% The tail durations will be t(end)-tailList
tailList = 1:69;

matFiles = dir([inDir '/*.mat']);

tic;
for i=1:length(matFiles)
    curFile = matFiles(i).name;
    load(fullfile(inDir,curFile));
    % 'Ct','Cp','Crr','t','maskCt','maskCrr','maskCp'
    Crr(Crr<0)=0;
    Ct(Ct<0)=0;
    Cp(Cp<0)=0;
    %%
    qtMask = max(Ct) > 0.01;
    numGoodVox = sum(qtMask(:));
    numVox = sum(maskCt(:));
    Ct = Ct(:,qtMask);
    maskCt(maskCt) = qtMask;
    %%
    [pkCE, ~, estKepRR, pkERRM] = CERRM(Ct,Crr,t);
    %% Loop over different starting frame for the tail
    for fTail=tailList
        estKtRR(i,fTail) = RRIFT(Cp(fTail:end),Crr(fTail:end),t(fTail:end),estKepRR);
        estKepRRS(i,fTail) = estKepRR;
    end
end
toc
%% Plot the results
% We will use the tail from frame 33 to frame 70 (end) as a target for comparison
refFrame = 33; 
% Plot the percent change in estimated KtransRR
y = PercentError(estKtRR,estKtRR(:,refFrame));
plot(t(end-1)-t(refFrame:end-1),y(:,refFrame:end))
ylim([-100 100])
ylabel('Percent Change in KtransRR')
xlabel('Duration of Ttail')
%% Find when error is above 20%
disp('The percent change exceed 20% when the tail duration is below (in minutes):')
disp(t(end) - t(find(any(abs(y)>20),1)));