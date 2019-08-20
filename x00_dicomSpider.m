%% Description
% This script renames the DICOM folders, e.g.
% before:   1.3.6.1.2.435.2345.453.45234253.4524353425234
% after:    3D DCE T1 MAPPING X 5
% The new TCIA client might be able to do the renaming, so this script is
% probably unnecessary. It's here for reference, in any case.

%% Initial setup
clearvars
addpath('./mfiles')

% Config: Set this to wherever the DICOM files are stored
startDir = '.\data\TCGA-GBM';

%% Pass the responsibility to another function
assert(isdir(startDir), 'Starting directory not found.')
SpiderFunc(startDir)
