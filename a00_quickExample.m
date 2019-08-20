% A quick single-voxel simulation showing how RRIFT works

clearvars

addpath('./mfiles');

%% Setup

rng(123);

% Set the temporal resolution
% Simulated data is initially generated at a higher temporal resolution...
initialTemporalRes = 0.1; % in seconds
% ... and then downsampled prior to fitting
temporalRes = 5; % in seconds

scanDuration = 10; % in minutes
bolusArrivalTime = 1; % in minutes

% Set the start/end time of the tail region
tailStartTime = scanDuration/2; % start time should be washout portion of AIF
tailEndTime = scanDuration;

% Set noise
noise = 0.02; % StdDev of noise in tissue of interest (units = mM)
noiseRR = noise/10; % StdDev of noise in reference tissue (units = mM)

% Define the tissue of interest parameters
ktrans = 0.25; % units: 1/min
ve = 0.4; % units: fractional
kep = ktrans/ve; % units: 1/min
vp = 0.05; % units: fractional

% Define the reference region parameters
ktransRR = 0.07; % units: 1/min
veRR = 0.14; % units: fractional
kepRR = ktransRR/veRR; % units: 1/min

%% Generate/Simulate Data

% Time should be in minutes
t = initialTemporalRes:initialTemporalRes:scanDuration*60;
t = t'/60;

% Generate the input function
Cp = GeorgiouAif(t,bolusArrivalTime);

% Simulate tissue curves
Ct = ToftsKety(Cp,[ktrans kep vp],t);
Crr = ToftsKety(Cp,[ktransRR kepRR],t);

% Add noise
Cp = Cp + noise * randn(size(Cp));
Ct = Ct + noise * randn(size(Ct));
Crr = Crr + noiseRR * randn(size(Crr));

% Downsample to desired temporal resolutions
downFactor = temporalRes/initialTemporalRes;
t = downsample(t, downFactor);
Cp = downsample(Cp, downFactor);
Ct = downsample(Ct, downFactor);
Crr = downsample(Crr, downFactor);

%% OVERVIEW
% We have four terms:
% - t: the time at each DCE frame, in minutes
% - Cp: the concentration in blood plasma, at each frame
% - Ct: the concentration in tissue of interest (tumour), at each frame
% - Crr: the concentration in the reference region (muscle), at each frame
% We will try to approaches: (i) Tofts model, (ii) Reference Region and Input Function Tail (RRIFT)

%% Tofts Model Fit

% Tofts model requires: Ct, Cp, t

estTofts = Tofts_LLSQ(Ct,Cp,t,1);
% (LLSQ is a poorly named function that performs linear fit of Tofts model)
% estTofts consists of: [ktrans, kep, vp]
% We will re-arrange this to be [ktrans, ve, vp] next:
estTofts = [estTofts(1), estTofts(1)/estTofts(2), estTofts(3)];

%% RRIFT

% RRIFT requires: Ct, Crr, t, and only the tail portion of Cp

% First fit the CERRM (requires Ct, Crr, t)
[estCERRM, ~, estKepRR] = CERRM(Ct, Crr, t);
% estCERRM consists of:
% [ktrans/ktransRR, ve/veRR, kep, vp/ktransRR, kepRR]
% We are only interested in the relative estimates, so we'll pick the
% first, second, and third element only
estCERRM = estCERRM([1,2,3]);
% i.e. estCERRM now only contains: [ktrans/ktransRR, ve/veRR, vp/ktransRR]

% Our next goal is to estimate ktransRR and veRR using only the tail portion of
% the input function. This is performed by:
% (Setup section defined start/end of tail)
tTail = t(t>tailStartTime & t<tailEndTime);
CpTail = Cp(t>tailStartTime & t<tailEndTime);
CrrTail = Crr(t>tailStartTime & t<tailEndTime);

% Below is implementation of Eq.8 in manuscript
num = CrrTail - CrrTail(1) + estKepRR * cumtrapz(tTail,CrrTail);
denum = cumtrapz(tTail,CpTail);
estKtransRR = denum\num;

% Estimating veRR is trivial since we have already estimated ktransRR and kepRR
estVeRR = estKtransRR/estKepRR;

% Now we convert relative estimates from CERRM to absolute estimates
estRRIFT = [estCERRM(1)*estKtransRR, estCERRM(2)*estVeRR, estCERRM(3)*estKtransRR];

%% Results

disp('===========================');

disp('Simulated ground truth [ktrans, ve, vp]')
disp([ktrans, ve, vp]);

disp('Tofts model estimates [ktrans, ve, vp]')
disp(estTofts)

disp('RRIFT estimates [ktrans, ve, vp]')
disp(estRRIFT)

disp('*******************');

disp('Simulated reference tissue [ktransRR, veRR]')
disp([ktransRR, veRR])

disp('RRIFT estimates [ktransRR, veRR]')
disp([estKtransRR, estVeRR])

disp('===========================');