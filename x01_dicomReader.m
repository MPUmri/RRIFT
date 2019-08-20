% Reads DICOM data for the Variable Flip Angle and DCE acquisitions, 
% and then does T1 mapping followed by signal-to-concentration conversion.
% The data is then saved for future processing.
% This function does not overwrite existing files, so delete any old files
% when re-running this script.

% The output of this script is available on OSF, so this script can be skipped

% Estimated runtime: ~8 minutes
% Output directory size: ~1.5 Gb
%% Initial setup
clearvars
addpath('./mfiles')

rootDir = './data/TCGA-GBM/';
% The filesystem at `rootDir` should resemble:
% rootDir
% +-- Patient 1
% |----+ Visit 1
% |    |----+ DCE
% |    |    |---- DCM1
% |    |    |---- DCM2
% |    |    |---- ...
% |    |----+ VFA
% |         |---- DCM1
% |         |---- DCM2
% |         |---- ...
% +-- Patient 2
% |----+ Visit 1
% |    |----+ DCE
% |    |    |---- DCM1
% |    |    |---- DCM2
% |    |    |---- ...
% |    |----+ VFA
% |         |---- DCM1
% |         |---- DCM2
% |         |---- ...

% The signal (DCE), T1 map (T1) and concentration (Ct) will be saved here:
outDir = './data/TCGA-GBM-Mat/';
outDirDCE = fullfile(outDir,'DCE');
outDirT1 = fullfile(outDir,'T1');
outDirCt = fullfile(outDir,'Ct');
outDirHDR = fullfile(outDir,'hdr');

if ~exist(outDir,'dir')
    mkdir(outDir)
    mkdir(outDirDCE)
    mkdir(outDirCt)
    mkdir(outDirT1)
    mkdir(outDirHDR)
end

%% Process the patient data
patientDir = dir(rootDir);
patientDir(1:2) = []; % Drop the first two entries because they are '.' and '..'

tic;
for p = 1:length(patientDir) 
    curPatient = patientDir(p).name % Display current patient in console
    patientPath = fullfile(rootDir,curPatient);
    
    % Each patient can have multiple visits
    visitDir = dir(patientPath);
    visitDir(1:2) = [];
    for q = 1:length(visitDir)
        curName = [curPatient '-' num2str(q) '.mat'];
        if exist(fullfile(outDirDCE,curName), 'file')
            % Skip if DCE output file already exists
            continue
        end
        visitPath = fullfile(patientPath, visitDir(q).name);
        dceDir = dir([visitPath '/*DYN*']);
        vfaDir = dir([visitPath '/*T1 MAP*']);
        if isempty(dceDir)
            error('Could not find DCE data')
        end
        dcePath = fullfile(visitPath, dceDir(1).name);
        vfaPath = fullfile(visitPath, vfaDir(1).name);
        [ctData, t1Data, m0Data, dceData, t, dceHdr, vfaHdr, flipAngles, flipData] = SortItOutForMe(dcePath, vfaPath);
        % Save as 'single' instead of 'double' to save disk space & faster loading
        ctData = single(ctData); t1Data = single(t1Data); 
        m0Data = single(m0Data); dceData = single(dceData);
        save(fullfile(outDirDCE,curName),'dceData','dceHdr','t');
        save(fullfile(outDirT1,curName),'t1Data','m0Data','vfaHdr','flipAngles','flipData');
        save(fullfile(outDirCt,curName),'ctData','t');
        save(fullfile(outDirHDR,curName),'dceHdr');
    end
end
toc
disp('Done reading DICOM data and saving results')
