%% Initialize
addpath('./mfiles')
clearvars

rng(12345)
load('./data/simMap.mat');
sigmaC = 0.03; % StdDev of noise (in units of mM)
TRes = 5; % Temporal resolution (in seconds)

% Properties for the simulated reference tissue
ktRR = 0.07; % Ktrans of ref.tissue, in units of 1/min
kepRR = 0.5; % kep of ref.tissue, in units of 1/min
veRR = ktRR/kepRR;
Crr = ToftsKety(Cp,[ktRR,kepRR],t);
%% Downsample and Add noise
Ct = downsample(simCt, TRes./initTRes);
Cp = downsample(Cp, TRes./initTRes);
Crr = downsample(Crr, TRes./initTRes);
t = downsample(t, TRes./initTRes);

[sT sX sY] = size(Ct);
Ct = reshape(Ct,[sT sX*sY]);
Ct = Ct + sigmaC * randn(size(Ct));
Cp = Cp + sigmaC * randn(size(Cp)) / (1-0.35);
Crr = Crr + 0.1 * sigmaC * randn(size(Crr));
%% Pick four voxels and fit them
% Chosen voxels (tried to avoid overlap of curves in figure)
idx = [5, 28, 53, 98];

% Fit the chosen voxels with CERRM
[pkCE, ~, estKepRR] = CERRM(Ct(:,idx),Crr,t);

% Use RRIFT, with AIF tail starting at 3 minutes into acquisition
fTail = find(t>3, 1);
[estKtRR, num, denum] = RRIFT(Cp(fTail:end), Crr(fTail:end), t(fTail:end), estKepRR);
%% Make figure

figure('Position',[300,300,1500,600])

% Muscle curve
subplot(2,2,1)
plot(t, Crr, 'g', 'LineWidth', 2)
ylim([-0.01 0.25])
xlabel('Time [min]');
ylabel('Concentration [mM]')
title('Concentration in Reference Region (Muscle)')
customizeFig(gca);

% Tissue curves
subplot(2,2,2)
plot(t, Ct(:,idx), 'LineWidth', 2)
ylim([-0.1 1.5])
xlabel('Time [min]');
ylabel('Concentration [mM]')
title('Concentration in Tissues')
customizeFig(gca);

% Blood plasma curve
subplot(2,2,3)
plot(t, Cp, 'r', 'LineWidth', 2)
ylim([-0.1 12])
xlabel('Time [min]');
ylabel('Concentration [mM]')
title('Concentration in Blood Plasma')
customizeFig(gca);

% RRIFT fit
subplot(2,2,4)
scatter(denum, num, 80,'r','LineWidth',2,'MarkerFaceColor',[1 1 1])
hold on;
plot(denum, denum*(denum\num), 'k-.', 'LineWidth', 2);
hold off;
ylim([-0.05 0.5])
ylabel('Numerator')
xlabel('Denominator')
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh, ['R^2: ' num2str(corr2(denum,num).^2)],'location','southeast')
legend boxoff
title('RRIFT Fit')
customizeFig(gca);
%% EXTRA: Fit all 100 combinations
% The above is very artificial because only four voxels were fitted. 
% This was done to make a simplified figure in the manuscript.
% In practice, tumour will have >100 voxels, and RRIFT should be more reliable. 
% This is done here where all 100 simulated voxels are fitted.
[~, ~, estKepRR_all] = CERRM(Ct,Crr,t);
fTail = find(t>3, 1);
[estKtRR_all, num, denum] = RRIFT(Cp(fTail:end), Crr(fTail:end), t(fTail:end), estKepRR_all);
%%
disp('')
disp('=========================')

disp('------Four-voxel fit-----')
disp('True KtransRR:')
disp(ktRR)
disp('Estimated KtransRR:')
disp(estKtRR)
disp('%Error in estimate')
disp(PercentError(estKtRR, ktRR))

disp('----All (100) voxel fit---')
disp('True KtransRR:')
disp(ktRR)
disp('Estimated KtransRR:')
disp(estKtRR_all)
disp('%Error in estimate')
disp(PercentError(estKtRR_all, ktRR))

%% MORE EXTRA: Try reference tissue plus vessel (RTPV) approach
% Note: re-run script with `sigmaC=0` to get values reported in Discussion
estVeRR_RRIFT = estKtRR./estKepRR;
estVeRR_RTPV = Crr(end)./Cp(end);
disp('')
disp('=========================')
disp('----Comparison to RTPV----')
disp('True veRR:')
disp(veRR)
disp('Estimated veRR from RRIFT (and % error):')
disp([estVeRR_RRIFT, PercentError(estVeRR_RRIFT, veRR)])
disp('Estimated veRR from RTPV (and % error):')
disp([estVeRR_RTPV, PercentError(estVeRR_RTPV, veRR)])

% Additional notes: The RTPV result improve if veRR is estimated by:
% >> mean(Crr(fTail:end)./Cp(fTail:end);
% or by:
% >> Cp(fTail:end)./Crr(fTail:end);
% particularly when noise is included (sigmaC =/= 0).
% For simplicity, we reported values for the noiseless case, in which case
% Crr(end)./Cp(end) should give an estimate for veRR according to the RTPV theory. 
%%