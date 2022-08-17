% Make Json copies of existing spreadsheet config files.
%
% Author(s): Rob Oldaker
%
% Modified by Yitong Li 
% The script is changed to a function, which can automatically find the
% file path.

function ConvertExcelFile2JsonFile(FileName)

%%
% Change the matlab path to the file path
RootPath = mfilename('fullpath');        % Get the path of this file
[RootPath,~,~]  = fileparts(RootPath);
[RootPath,~,~]  = fileparts(RootPath);
[RootPath,~,~]  = fileparts(RootPath);
cd(RootPath);                            % Change the current address

%%
FilePath = fileparts(which(FileName));

% Check if the file name is proper
if isempty(FilePath)
    error(['User excel could not be found. Please double check if the file name and file type are correct (".xls", ".xlsm", ".xlsx", etc) or if the file path is added.'])
end

% Get the folder name of the excel file
FolderName = erase(FilePath,RootPath);

% Get the full name of the excel file
FileFullName = [FolderName,'\',FileName];

% Remove '\' at the begining of the file name if needed.
if strcmp(FileFullName(1),'\')
    FileFullName = FileFullName(2:end);
end

%%
% Convert excel file to json file
SimplusGT.Toolbox.Excel2Json(FileFullName);

end