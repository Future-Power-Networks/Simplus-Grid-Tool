% Please read the comments in this file carefully, and use this file to run
% toolbox.

%% Tips
%
% Please read manuals in the "Documentations" folder if you want to know
% more details about this tool.
%
% Please ensure that the toolbox is installed first, by running
% "InstallSimplusGT.m" once.
%
% For changing default user data, please use "UserData.xlsx". More examples
% can be found in "Examples" folder.

%% Clear matlab
clear all;  % Clear Matlab workspace
clc;        % Clear Matlab command window
close all;  % Close all figures, etc

%% Set user data
% "UserData.xlsx" and "UserData.json" contain the data of an example 4-bus
% generator-inverter-composite power system. Please feel free to change
% them. 
% 
% ".xlsm" or ".xlsx" is the excel file, and ".json" is the corresponding json
% file. Users can easily convert an Excel file to a json file by calling
% this function: 
SimplusGT.Toolbox.ConvertExcelFile2JsonFile('UserData.xlsm');
UserData = 'UserData.json';

%%
% Other example power systems (in "Examples" folder):
%
% Pure ac power system examples:
% UserData = 'SgInfiniteBus.json';              % Single synchronous generator and infinite bus
% UserData = 'GflInverterInfiniteBus.json';   	% Single grid-following inverter and infinite bus
% UserData = 'GfmInverterInfiniteBus.json';   	% Single grid-forming inverter and infinite bus
% UserData = 'IEEE_14Bus.json';
% UserData = 'IEEE_30Bus.json';
% UserData = 'IEEE_57Bus.json';
% UserData = 'NETS_NYPS_68Bus.json';
%
% Pure dc power system examples:
% UserData = 'GfdBuckInfiniteBus.json';         % Single grid-feeding buck converter and infinite bus
%
% Hybrid ac-dc power system examples:
% UserData = 'Hybrid_test_v1.json';             % A 4-bus hybrid ac-dc system

%% Run toolbox
InputData = SimplusGT.JsonDecoder(UserData);
SimplusGT.Toolbox.Main();