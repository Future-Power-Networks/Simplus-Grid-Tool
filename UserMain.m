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
% UserDataName = 'BessInfiniteBus';             % Single battery energy storage system and infinite bus
% UserDataName = 'PVGfmInfiniteBus';            % Single Photovoltaic(GFM) and infinite bus
UserDataName = 'PVGflInfiniteBus';            % Single Photovoltaic(GFL) and infinite bus
% UserDataName = 'IEEE_14Bus';
% UserDataName = 'IEEE_30Bus';
% UserDataName = 'IEEE_57Bus';
% UserDataName = 'AU14Gen_59Bus';
% UserDataName = 'NETS_NYPS_68Bus';
% UserDataName = 'BESS_Plant_10Bus';            % A 10-unit battery energy storage plant
% UserDataName = 'PV_Plant_10Bus';              % A 10-unit photovoltaic plant
% UserDataName = 'PV_BESS_Hybrid_Plant_10Bus';  % A 10-unit photovoltaic and energy storage hybrid plant
%
% Dc power system examples:
% UserDataName = 'GfdBuckInfiniteBus';         % Single grid-feeding buck converter and infinite bus
% UserDataName = 'TwoBusGfdBuck';              % Two buck converters
%
% Hybrid ac-dc power system examples:
% UserDataName = 'Hybrid_4Bus';             % A 4-bus hybrid ac-dc system
% UserDataName = 'Hybrid_28Bus';
% UserDataName = 'HVDC_Infbus_4Bus';        % HVDC system connected to inf buses
% UserDataName = 'HVDC_SG_4Bus';            % HVDC system connected to equivalent SG buses
% UserDataName = 'MTDC_Infbus_4Bus';        % MTDC system connected to inf buses

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