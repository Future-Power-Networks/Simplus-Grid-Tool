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
% file = 'UserData.xlsm';
% SimplusGT.Toolbox.Excel2Json(file);
% 
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
% 
% file = 'Examples\AcPowerSystem\SingleApparatusInfiniteBus\GflInverterInfiniteBus.xlsx';
% SimplusGT.Toolbox.Excel2Json(file);
% 
% file = 'Examples\AcPowerSystem\SingleApparatusInfiniteBus\GfmInverterInfiniteBus.xlsx';
% SimplusGT.Toolbox.Excel2Json(file);
% 
file = 'Examples\AcPowerSystem\NET_NYPS_68Bus\NETS_NYPS_68Bus.xlsx';
SimplusGT.Toolbox.Excel2Json(file);

file = 'Examples\NETS_NYPS_68Bus_LowDamping.xlsx';
SimplusGT.Toolbox.Excel2Json(file);

file = 'Examples\NETS_NYPS_68Bus_HighDamping.xlsx';
SimplusGT.Toolbox.Excel2Json(file);

file = 'Examples\NETS_NYPS_68Bus_LowDamping_SingleSgCase.xlsx';
SimplusGT.Toolbox.Excel2Json(file);

file = 'Examples\NETS_NYPS_68Bus_HighDamping_SingleSgCase.xlsx';
SimplusGT.Toolbox.Excel2Json(file);
% 
% file = 'Examples\AcPowerSystem\IEEE_57Bus\IEEE_57Bus.xlsx';
% SimplusGT.Toolbox.Excel2Json(file);
% 
% file = 'Examples\AcPowerSystem\IEEE_30Bus\IEEE_30Bus.xlsx';
% SimplusGT.Toolbox.Excel2Json(file);
% 
% file = 'Examples\AcPowerSystem\IEEE_14Bus\IEEE_14Bus.xlsx';
% SimplusGT.Toolbox.Excel2Json(file);

% file = 'Examples\TestSynchronisation\Test_68Bus_IBR.xlsx';
% SimplusGT.Toolbox.Excel2Json(file);
% file = 'Examples\TestSynchronisation\Test_68Bus_IBR_17.xlsx';
% SimplusGT.Toolbox.Excel2Json(file);
% file = 'Examples\TestSynchronisation\Test_68Bus_IBR_17_14.xlsx';
% SimplusGT.Toolbox.Excel2Json(file);