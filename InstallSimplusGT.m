% Read me please:
% This file installs the toolbox on the user PC. Users do not need to read
% or understand the codes in this file. Please just run it directly and run
% it only ONCE. Please use "UserMain.m" in the root path to run the toolbox
% after installation.

% Author(s): Yitong Li, Yunjie Gu

%%
% Notes for developers:
%
% The simulink library may need to be updated programmatically here, when
% adding Simplus Library.

%%
clear all
clc
close all

%%
% Change the current folder
fprintf('Changing the current folder to the toolbox folder...\n')
MfilePath = mfilename('fullpath');
[RootPath,~,~]  = fileparts(MfilePath);
cd(RootPath);

% Check the matlab version
fprintf('Checking the Matlab version...\n')
MatlabVersion = version('-release');
MatlabVersion = MatlabVersion(1:(end-1));
MatlabVersionYear = MatlabVersion;
if str2double(MatlabVersionYear)<2015
    error(['Error: Please use Matlab version 2015a or later!']);
end

% Check if a previous version of SimplusGT has been installed
fprintf('Checking if the toolbox has been installed before...\n')
if exist('SimplusGT')~=0
    error(['Error: The toolbox has been installed on this PC/laptop before. Please unstall the old version first!']);
end

% Add folder to path
fprintf('Installing...\n')
addpath(genpath([RootPath,'/Examples']));         	% Add "Examples" folder
addpath(genpath([RootPath,'/Library']));          	% Add "Library" folder
addpath(genpath([RootPath,'/Documentations']));   	% Add "Documentations" folder
addpath(genpath(RootPath));                         % Add root path
savepath;

% Convert the toolbox lib to the required version
fprintf('Converting the toolbox library to the required Matlab version, please wait a second...\n')
warning('off','all')    % Turn off all warnings
open_system('SimplusGT_2015a.slx');
save_system('SimplusGT_2015a.slx','Library/SimplusGT.slx');
close_system('SimplusGT.slx');
warning('on','all')     % Turn on all warnings
clc

%%
% % Installation is completed
% DlgTitle = 'Congratulations!';
% DlgQuestion = 'The toolbox is installed successfully! Do you want to run "UserMain.m" to experience it NOW?';
% choice = questdlg(DlgQuestion,DlgTitle,'Yes','No','Yes');
% if strcmp(choice,'Yes')
%     open('UserMain.m');
%     run('UserMain.m');
% else
%  	msgbox('The installation is completed. Please run "UserMain.m" later in the root path for using the toolbox.');
% end

%%
fprintf('Congratulations!\n')
fprintf('The installation is completed.\n')
fprintf('Please run "UserMain.m" later in the root path for using the toolbox!\n')
fprintf('Please read manuals in the "Documentations" folder if you want to know more details.\n')