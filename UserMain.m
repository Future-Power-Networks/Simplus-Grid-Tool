%% Readme
%
% Please ensure that the toolbox is installed first by
% "InstallSimplusGT.m".
%
% Default user data is saved in "UserData.xlsm" and "UserData.json". More
% examples can be found in "Examples" folder.
%
% Please read manuals in the "Documentations" folder if you want to know
% more details about this tool.

%% Clear matlab
clear all;  % Clear Matlab workspace
clc;        % Clear Matlab command window
% close all;  % Close all figures, etc

%% Change the current folder of matlab
cd(fileparts(mfilename('fullpath')));

%% Set user data name
UserDataType = 1;
% UserDataName = 'Data14Bus';
% UserDataName = 'Data68Bus';
% UserDataName = 'Data68Bus_VirtualRC';
UserDataName = 'Data118Bus';
% UserDataName = 'Data118Bus_VirtualRC';
% UserDataName = 'Data118Bus_Sg';

%% Run toolbox
SimplusGT.Toolbox.Main();

%% Find mode
% ModeIndex = intersect(find(imag(EigVecHz)<40),find(imag(EigVecHz)>20)) 
% EigVecHz(ModeIndex)