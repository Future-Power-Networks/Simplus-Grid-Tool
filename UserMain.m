%% Readme
%
% Default 4-bus power system user data is saved in "UserData.xlsm" and
% "UserData.json". More examples can be found in "Examples" folder.
% (默认用户数据见"UserData.xlsm"和"UserData.json".更多案例系统见
% "Examples"文件夹.)
%
% More manuals are available in the "Documentations" folder.
% (更多手册见"Documentations"文件夹.)

%% Clear matlab
clear all; clc; close all; 

%% Other examples
UserDataName = 'UserData';      % Default 4-bus system

% Other example power systems in "Examples" folder:
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

%% Change the current folder of matlab
cd(fileparts(mfilename('fullpath')));

%% Set user data type
% If user data is in excel format, please set 1. If it is in json format,
% please set 0.
%(若用户数据是Excel形式,请设置1.若用户数据为json形式，请设置0。）
UserDataType = 1;

%% Run toolbox
SimplusGT.Toolbox.Main();