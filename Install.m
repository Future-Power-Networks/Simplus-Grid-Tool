% Author(s): Yitong Li

%%
clear all
clc
close all

%%
% Change the current folder
fprintf('Changing the current folder to the toolbox folder...\n')
MfilePath = mfilename('fullpath');
[RootPath,~,~]  = fileparts(MfilePath);
cd(RootPath);

% Add folder to path
fprintf('Installing...\n')
addpath(genpath([RootPath,'/Library']));
addpath(genpath([RootPath,'/MultiInverterCase_Simulation']));
% addpath(genpath(RootPath));                         % Add root path
savepath;

% Print
fprintf('The toolbox is installed! \n')