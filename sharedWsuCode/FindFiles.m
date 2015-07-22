function fileCell = FindFiles(pathStr)
%returns files in a cell array found using path Str - which may contain wildcards


%unix systems only - only wildcard accepted is *?
% fileMat=ls(pathStr); %returns the paths in a matrix, so then we must separate and remove white spaces
% 
% fileCell=textscan(fileMat, '%s', 'multipleDelimsAsOne', true);
% fileCell= fileCell{1};

%decompose pathStr
[folderStr fileStr extStr]=fileparts(pathStr);
if isempty(strfind(folderStr, '*'))
%there are no wildcards in the folder path, we can just do a dir
    fileCell=dir(pathStr);
    fileCell={fileCell.name};
end

folderList={};
while ~isempty(strfind(folderStr, '*'))    
    [folderStr thisFolder]=fileparts(folderStr);
    folderList=[thisFolder folderList];
end
%now go through and collect filenames starting at earliest wildcard

foldersToSearch=[];
while ~isempty(folderList)
%     folderStr=[folderStr filesep folderList{1}];
    foundFolders=dir([folderStr filesep folderList{1}]);
    foundFolders={foundFolders.name};
    folderList{1}=[];
    for fIndx=1:length(foundFolders)
        foldersToSearch{end+1}=[folderStr filesep foundFolders{fIndx}];        
    end
    
end