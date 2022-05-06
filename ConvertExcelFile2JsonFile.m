% Make Json copies of existing spreadsheet config files.
%
% Author(s): Rob Oldaker, Yitong Li

%%
% Clear
clear all
clc
close all

%%
% Change the matlab path to the file path
PathStr = mfilename('fullpath');        % Get the path of this file
[PathStr,~,~]  = fileparts(PathStr);
cd(PathStr);                            % Change the current address

%%
% 
file = 'Examples\AcPowerSystem\NET_NYPS_68Bus\NETS_NYPS_68Bus.xlsx';
SimplusGT.Toolbox.Excel2Json(file);

file = 'Examples\68Bus_HighInertia_InterAreaMode.xlsx';
SimplusGT.Toolbox.Excel2Json(file);

file = 'Examples\68Bus_MedInertia_InterAreaMode.xlsx';
SimplusGT.Toolbox.Excel2Json(file);

file = 'Examples\68Bus_LowInertia_InterAreaMode.xlsx';
SimplusGT.Toolbox.Excel2Json(file);

file = 'Examples\68Bus_HighInertia_LocalMode.xlsx';
SimplusGT.Toolbox.Excel2Json(file);

file = 'Examples\68Bus_MedInertia_LocalMode.xlsx';
SimplusGT.Toolbox.Excel2Json(file);

file = 'Examples\68Bus_LowInertia_LocalMode.xlsx';
SimplusGT.Toolbox.Excel2Json(file);