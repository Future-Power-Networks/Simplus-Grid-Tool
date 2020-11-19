% This function installs the toolbox into the customer PC.

%% Notes
%
% After running this file, simply run "CustomerMain.m" to run the toolbox.

%%
clear all
clc
close all

%% Add folder to path
addpath(genpath(pwd));  % Add path
savepath;               % Save path

%%
fprintf('Install SimplexPowerSystem successfully! \n')
fprintf('Please run "CustomerMain.m" to run the toolbox.')