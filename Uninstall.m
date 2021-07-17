clear all
clc
close all

% Change the current folder
mfile_name = mfilename('fullpath');
[pathstr,~,~]  = fileparts(mfile_name);
cd(pathstr);

% Remove folder from path
rmpath(genpath(pwd));
restoredefaultpath;
clc

fprintf('The toolbox is uninstalled sucessfully! \n')