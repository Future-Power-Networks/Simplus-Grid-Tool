% Author(s): Yitong Li

%%
clear all
clc
close all

%% 
% Change the current folder
fprintf('Changing the current folder to the toolbox folder...\n')
mfile_name = mfilename('fullpath');
[pathstr,~,~]  = fileparts(mfile_name);
cd(pathstr);

% Remove folder from path
rmpath(genpath(pwd));
restoredefaultpath;
savepath;
clc;

% Print
fprintf('The toolbox is uninstalled! \n')