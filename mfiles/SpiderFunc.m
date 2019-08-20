function SpiderFunc(startDir)
% Crawls through directories looking for dicom files
% If it finds a dicom file, it renames the folder to series description.

% Alternative description: Does whatever a SpiderFunc does.

% Search for dicom files in given directory
if ~isempty(dir(fullfile(startDir, '*.dcm')))
    dcmFiles = dir(fullfile(startDir, '*.dcm'));
    % Read first file to get SeriesDescription
    dcmFile = dcmFiles(1).name;
    dcmInfo = dicominfo(fullfile(startDir,dcmFile));
    sDescription = dcmInfo.SeriesDescription;
    % Remove special characters which can lead to bad folder name
    sDescription = regexprep(sDescription,'[*:/?.\\<>]',' ');
    [pathstr,~,~] = fileparts(startDir); % Get current path
    newDir = fullfile(pathstr,sDescription); % Make new dir name
    if ~isempty(findstr(startDir, newDir)) % Why didn't I use strcmp here [?]
        % This will avoid renaming files that don't need it
        return
    end
    % Console will show progress
    startDir
    newDir
    movefile(startDir,newDir); % Rename folder
    return
end

% If no dicom file is found...

% Get list of folders in directory
fileList = dir(startDir);

% Search through all these folders for dicom
for i=3:length(fileList) % Ignore first two because they are '.' and '..'
    SpiderFunc(fullfile(startDir,fileList(i).name));
end
    
end

