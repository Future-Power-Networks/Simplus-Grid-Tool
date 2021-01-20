% This function installs the SimplexPowerSystem toolbox into the customer
% PC.

%% Notes

% Please run this file just ONCE.

%%
clear all
clc
close all

fprintf('Installing SimplexPowerSystem...\n')

%% Add folder to path
addpath(genpath(pwd));  % Add path
savepath;               % Save path

%% Save lib to the corresponding matlab version
%MatlabVersion = ver('MATLAB').Release;
MatlabVersion = version('-release');
MatlabVersion = MatlabVersion(1:(end-1));
MatlabVersionYear = MatlabVersion;
if str2double(MatlabVersionYear)<2015
    error(['Error: Please use Matlab version 2015a or later!']);
end

fprintf('Converting the toolbox library to the required Matlab version, please wait a second...\n')
warning('off','all')    % Turn off all warnings
open_system('SimplexPS_2015a.slx');
save_system('SimplexPS_2015a.slx','Library\SimplexPS.slx');
close_system('SimplexPS.slx');
warning('on','all')     % Turn on all warnings
clc

%%
% fprintf('SimplexPowerSystem is installed successfully! \n')
% fprintf('Do you want to run "CustomerMain.m" now to start using the toolbox?\n')

DlgTitle    = 'Congratulations!';
DlgQuestion = 'SimplexPowerSystem is installed successfully! Do you want to run "UserMain.m" now to experience the toolbox?';
choice = questdlg(DlgQuestion,DlgTitle,'Yes','No','Yes');

if strcmp(choice,'Yes')
   open('CustomerMain.m');
   run('CustomerMain.m');
else
    msgbox('The installation of SimplexPowerSystem is completed. Please run "UserMain.m" later in the root path when you want to experience the toolbox!');
end