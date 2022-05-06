% Author(s): Yitong Li

%% Clear matlab
clear all;  % Clear Matlab workspace
clc;        % Clear Matlab command window
close all;  % Close all figures, etc

%% Get input data
ConvertExcelFile2JsonFile();

% InputData = SimplusGT.JsonDecoder('NETS_NYPS_68Bus.json');

% InputData = SimplusGT.JsonDecoder('68Bus_HighInertia_InterAreaMode.json');
InputData = SimplusGT.JsonDecoder('68Bus_MedInertia_InterAreaMode.json');
% InputData = SimplusGT.JsonDecoder('68Bus_LowInertia_InterAreaMode.json');

% InputData = SimplusGT.JsonDecoder('68Bus_HighInertia_LocalMode.json');
% InputData = SimplusGT.JsonDecoder('68Bus_MedInertia_LocalMode.json');
% InputData = SimplusGT.JsonDecoder('68Bus_LowInertia_LocalMode.json');

%% Run toolbox
SimplusGT.Toolbox.Main();

%% Set C parameters of the LC filter
Rcap = 0.0001;
Ccap = 0.03/Wbase;