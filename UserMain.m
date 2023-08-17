%% Readme
%
% Default 4-bus power system user data is saved in "UserData.xlsm" and
% "UserData.json". More examples can be found in "Examples" folder.
%
% More manuals are available in the "Documentations" folder.

%% Clear matlab
clear all; clc; close all; 

%% User data
UserDataName = 'UserData';      % Default 4-bus system

% Example power systems in "Examples" folder:
%
% Ac power system examples:
% UserDataName = 'SgInfiniteBus';               % Single synchronous generator and infinite bus
% UserDataName = 'GflInverterInfiniteBus';   	% Single grid-following inverter and infinite bus
% UserDataName = 'GfmInverterInfiniteBus';   	% Single grid-forming inverter and infinite bus
% UserDataName = 'IEEE_14Bus';
% UserDataName = 'IEEE_30Bus';
% UserDataName = 'IEEE_39Bus';
% UserDataName = 'IEEE_57Bus';
% UserDataName = 'AU14Gen_59Bus';
% UserDataName = 'NETS_NYPS_68Bus';
%
% Dc power system examples:
% UserDataName = 'GfdBuckInfiniteBus';         % Single grid-feeding buck converter and infinite bus
%
% Hybrid ac-dc power system examples:
% UserDataName = 'Hybrid_test_v1';             % A 4-bus hybrid ac-dc system

%% Change the current folder of matlab
cd(fileparts(mfilename('fullpath')));

%% Set user data type
% If user data is in excel format, please set 1. If it is in json format,
% please set 0.
UserDataType = 1;

%% Run toolbox
SimplusGT.Toolbox.Main();  

%% Matlab app
if 1; ModalAnalysisApp; end      % Modal analysis