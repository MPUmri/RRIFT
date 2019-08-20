%% Description
% This script does a few pre-processing steps, including:
% - Loads concentration data
% - Applies masks for tumour, muscle, and blood vessel
% - Collects the concentration-time curve for tumour/muscle/AIF
% - Calculates T1 for tumour, muscle, input function
% - Computes SNR for tumour, muscle, input function
% - Saves results as a .mat file for further analysis

% This script is essentially collecting the raw data without analyzing it.
% Estimated runtime: ~30 seconds

%% Initialize
clearvars;
fclose('all');
addpath('./mfiles');

%% Configuration

% Input directory containing the 'Ct' and 'DCE' directory
rootDir = './data/TCGA-GBM-Mat';

% Directories containing the masks
maskDirCp = './data/TCGA-GBM-Masks/AIF'; % Arterial input function
maskDirCrr = './data/TCGA-GBM-Masks/Muscle'; % Reference region
maskDirCt = './data/TCGA-GBM-Masks/Tumour'; % Tissue of interest

% Output directory for the formatted data
outDir = './data/TCGA-GBM-Results/c01_preprocessed';

% Overwrite existing files?
doOverwrite = true;

% Assumed value for blood hematocrit
hct = 0.4;

%% Main

% Create output directory if it doesn't exist
if ~exist(outDir,'dir')
    mkdir(outDir)
end

matFiles = dir([maskDirCt '/*.mat']);
tic;
for i=1:length(matFiles)
    curFile = matFiles(i).name;
    outFile = fullfile(outDir,curFile);
    if exist(outFile,'file') && ~doOverwrite
        continue
    end
    % Display the current patient name in console as a way of tracking progress
    curFile
    %% Load masks
    % Each mask .mat contains a variable named 'mask', so have to rename it.
    % The masks were saved as 'double' so also have to convert to 'logical'.
    load(fullfile(maskDirCp,curFile));
    maskCp = logical(mask);
    
    load(fullfile(maskDirCrr,curFile));
    maskCrr = logical(mask);
    
    load(fullfile(maskDirCt,curFile));
    maskCt = logical(mask);
    
    clearvars mask
    %% Load concentration data
    load(fullfile(rootDir,'Ct',curFile));
    % Provides:
    % - ctData: [4D array - X,Y,Z,time] concentration data
    % - t : [1D array - time] time in units of minutes
    
    % Melt 4D array to 2D array with dims: [time, X*Y*Z]
    % This makes it easier to process the data
    Ct = unravel(ctData)';
    
    clearvars ctData
    %% Get arterial input function, muscle curve, and tumour curves
    % AIF
    Cp = Ct(:,maskCp(:))./(1-hct);
    indCp = find(maskCp(:)>0,1); % This is used later to get T1 of arterial voxel
    % Reference region
    Crr = double(mean(Ct(:,maskCrr),2)); % Mean curve from muscle mask
    % Tissue of interest
    Ct = double(Ct(:,maskCt)); % All voxels in tumour mask
    %% Get T1 values
    load(fullfile(rootDir,'T1',curFile));
    % Provides:
    % vfaHDR : DICOM header for variable flip angle images
    % flipAngles : Acquired flip angles
    % flipData : Variable flip angle data (i.e. T1 weighted signal)
    % m0Data : Estimated magnetization (coefficient in SPGR equation)
    % t1Data : T1 map for all slices
    
    T1Cp = t1Data(indCp);
    T1Crr = mean(t1Data(maskCrr));
    T1Ct = mean(t1Data(maskCt));
    
    clearvars vfaHDR flipAngles flipData m0Data t1Data
    %% Get SNR and CNR
    load(fullfile(rootDir,'DCE',curFile));
    % Provides:
    % - dceData: 4D array containing MRI signal at each frame and slice
    % - dceHdr: DICOM header of the DCE acquisition
    % - t: time, in minutes, at each DCE frame
    
    % Melt 4D array to 2D array with dims: [time, X*Y*Z] to simplify next steps
    dceData = unravel(dceData)';
    
    % Calculate global SNR
    preContrastSignal = double(dceData(1:5,maskCt));
    snr = nanmean( nanmean(preContrastSignal) ./ nanstd(preContrastSignal));
    
    % Calculate global CNR
    bolusArrivalFrame = 5; % This might fail for one patient which has arrival at ~3rd frame
    preContrastCt = Ct(1:bolusArrivalFrame,:);
    % Remove extreme outliers (first and last percentile) so they don't corrupt stdDev
    preCtQt = quantile(preContrastCt(:),[.01 .99]);
    preContrastCt = preContrastCt(preContrastCt(:)>preCtQt(1) & preContrastCt(:)<preCtQt(2));
    sigmaCt = nanstd(preContrastCt(:));
    cnr = max(Cp)./sigmaCt; % CNR is defined as peak of input function / std deviation of concentration  
    %% Save results for further processing
    save(outFile,'Ct','Cp','Crr','t',...
    'maskCt','maskCrr','maskCp','T1Cp','T1Crr','T1Ct',...
    'cnr','snr','sigmaCt')
end
toc
disp('Done preprocessing data')