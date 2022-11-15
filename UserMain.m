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
close all;  % Close all figures, etc

%% Change the current folder of matlab
cd(fileparts(mfilename('fullpath')));

%% Set user data name
% "UserData.xlsm" and "UserData.json" contain the data of an example 4-bus
% generator-inverter-composite power system in excel and josn format
% respectively. Please feel free to change them.
% UserDataName = 'UserData';

%% Set user data type
% If user data is in excel format, please set 1. If it is in json format,
% please set 0.
UserDataType = 1;

%% Other examples
% Other example power systems (in "Examples" folder):
%
% Pure ac power system examples:
% UserDataName = 'SgInfiniteBus';               % Single synchronous generator and infinite bus
% UserDataName = 'GflInverterInfiniteBus';   	% Single grid-following inverter and infinite bus
% UserDataName = 'GfmInverterInfiniteBus';   	% Single grid-forming inverter and infinite bus
% UserDataName = 'IEEE_14Bus';
% UserDataName = 'IEEE_30Bus';
% UserDataName = 'IEEE_57Bus';
% UserDataName = 'NETS_NYPS_68Bus';
%
% Pure dc power system examples:
% UserDataName = 'GfdBuckInfiniteBus';         % Single grid-feeding buck converter and infinite bus
%
% Hybrid ac-dc power system examples:
% UserDataName = 'Hybrid_test_v1';             % A 4-bus hybrid ac-dc system

UserDataName = 'GfmInverterInfiniteBusTest';

%% Run toolbox
SimplusGT.Toolbox.Main();