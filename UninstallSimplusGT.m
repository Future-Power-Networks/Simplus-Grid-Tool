% This file uninstalls the SimplusGT.
% (运行此文件卸载SimplusGT)

% Author(s): Yitong Li

function UninstallSimplusGT()

clc

%% Generate bar msg
UninstallMsg = 'Uninstalling Simplus Grid Tool. This may take a few seconds...               ';
bar = waitbar(0.3,UninstallMsg,'Name','Simplus');
bar.CloseRequestFcn = '';

%% Uninstall
% Change the current folder
fprintf('Changing the current folder to the toolbox folder...\n')
MfileName = mfilename('fullpath');
[RootPath,~,~]  = fileparts(MfileName);
cd(RootPath);

% Remove folder from path
rmpath(genpath(pwd));
restoredefaultpath;
savepath;
clc;

fprintf('Simplus Grid Tool successfully uninstalled! \n')
fprintf('Many thanks for using! \n')
fprintf('Please delete the files manually if you want to fully remove Simplus Grid Tool from your PC. \n')

%% Remove bar msg
beep;
waitbar(1,bar,UninstallMsg,'Name','Simplus');
pause(1);
delete(bar);

msgbox('Simplus Grid Tool successfully unstalled!','Simplus');

end