% the function will save a model in 2015a format under 'install' folder
%
% Author(s): Yitong Li, Yunjie Gu
%
% parameter: 
% model = path/to/file/model.ext, in relative path to the
% current (working) path, and .ext is the extention name of the model,
% which is .slx by default if left empty
%
% usage: 
% method(1)
% cd('path/to/file')
% simplus.archive('model')
% method(2)
% simplus.archive('path/to/file/model')
% method(3)
% simplus.archive('path/to/file/model.mdl')

function Archive(model)

    current_path = pwd;

    [file_path,file_name,file_ext] = fileparts(model);
    if isempty(file_name)
        return;
    end
    if isempty(file_ext)
        file_ext = '.slx';
    end

    if ~isempty(file_path)
        cd(file_path);
    end
    if ~isfile([file_name file_ext])
        cd(current_path);
        return;
    end

    if ~isfolder('install')
        mkdir('install');
    end
    
    archive_path = fullfile(current_path,file_path,'Install',[file_name '_Archive' file_ext]); 
    load_system(file_name);
    save_system(file_name,archive_path,'ExportToVersion','R2015A');
    close_system(file_name);
    cd(current_path);
end

