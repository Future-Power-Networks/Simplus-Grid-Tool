%% Readme
%
% Default 4-bus power system user data is saved in "UserData.xlsm" and
% "UserData.json". More examples can be found in "Examples" folder.
%
% More manuals are available in the "Documentations" folder.

%% Clear matlab
clear all; clc; close all; 

%% User data
% UserDataName = 'UserData';      % Default 4-bus system

% Example power systems in "Examples" folder:
%
UserDataName = 'IEEE68BusSG';       % IEEE 68 bus system with synchronous machines
% UserDataName = 'IEEE68BusGFMGFL';   % IEEE 68 bus system with grid-forming and grid-following inverters


%% Change the current folder of matlab
cd(fileparts(mfilename('fullpath')));

%% Set user data type
% If user data is in excel format, please set 1. If it is in json format,
% please set 0.
UserDataType = 1;

%% Run toolbox
SimplusGT.Toolbox.Main();  

%% Matlab app
if 0; ModalAnalysisAPP; end      % Modal analysis