%% Readme

%% Clear matlab
clear all;  % Clear Matlab workspace
clc;        % Clear Matlab command window
% close all;  % Close all figures, etc

%% Change the current folder of matlab
cd(fileparts(mfilename('fullpath')));
UserDataType = 1;

%%
UserDataName = 'TestSingleGfl';
% UserDataName = 'TestSingleGfm';

%% Run toolbox
SimplusGT.Toolbox.Main();

