% To users:
% Please use this file to run toolbox.

%% Tips
%
% Please ensure that the toolbox is installed first, by running
% "InstallSimplexPS.m" once.
%
% The toolbox defaultly saves the results into Workspace, prints the key
% results in Command Window, and plots key figures.
%
% For changing default user data, please use "UserData.xlsx". More examples
% can be found in "Examples" folder.

%% Clear matlab
clear all;  % Clear matlab workspace
clc;        % Clear matlab command window
close all;  % Close all figures, etc

%% Set user data
UserData = 'Duality_14Bus.xlsx';

%% Run toolbox
SimplexPS.Toolbox.Main();