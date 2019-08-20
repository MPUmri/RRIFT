function [foundDirs] = findDir(startDir, targetString, foundDirs)
% This functions walks through a starting directory and searchs for folders
% with the named matching the 'targetString'
% Returns cell array containg path to found directories.
% Example usage: foundDirs = findDir('./', 'findThisString')

	% User should only provide two inputs, so we start with no found dirs    
    if nargin<3
        foundDirs = {};
    end

    % See if any folders in current directory match the target string
    dirList = dir(startDir);
    targetList = find(contains({dirList.name},targetString));
    if length(targetList)
        foundDirs{end+1}=startDir;
        return
    end
    
    % Then cycle through each subfolder to find matches
    for i=3:length(dirList)
        foundDirs = findDir(fullfile(startDir,dirList(i).name), targetString, foundDirs);
    end
end