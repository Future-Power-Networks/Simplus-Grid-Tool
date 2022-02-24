% This file uninstalls the SimplusGridTool toolbox.

% Author(s): Yitong Li

%%
clear all
clc
close all

%% 
% Change the current folder
fprintf('Changing the current folder to the toolbox folder...\n')
MfileName = mfilename('fullpath');
[RootPath,~,~]  = fileparts(MfileName);
cd(RootPath);

% Remove folder from path
rmpath(genpath(pwd));
restoredefaultpath;
savepath;
clc;

fprintf('The toolbox is uninstalled! \n')
fprintf('Many thanks for using! \n')
fprintf('Please delete the files manually if you want to fully remove this toolbox from your PC. \n')