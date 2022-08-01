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
% Convert excel file to json file
file = 'UserData.xlsm';
SimplusGT.Toolbox.Excel2Json(file);

% file = 'Examples\HybridPowerSystem\Hybrid_test_v1.xlsx';
% SimplusGT.Toolbox.Excel2Json(file);
% 
% file = 'Examples\HybridPowerSystem\Hybrid_test_v2.xlsx';
% SimplusGT.Toolbox.Excel2Json(file);
% 
% file = 'Examples\DcPowerSystem\SingleMachineInfiniteBus\GfdBuckInfiniteBus.xlsx';
% SimplusGT.Toolbox.Excel2Json(file);
% 
% file = 'Examples\AcPowerSystem\SingleApparatusInfiniteBus\SgInfiniteBus.xlsx';
% SimplusGT.Toolbox.Excel2Json(file);

file = 'Examples\AcPowerSystem\SingleApparatusInfiniteBus\GflInverterInfiniteBus.xlsx';
SimplusGT.Toolbox.Excel2Json(file);

file = 'Examples\AcPowerSystem\SingleApparatusInfiniteBus\GfmInverterInfiniteBus.xlsx';
SimplusGT.Toolbox.Excel2Json(file);

