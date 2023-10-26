% This file installs SimplusGT. Users do not need to read or understand the
% codes in this file. Please just run it directly and run it ONCE. Please
% use "UserMain.m" in the root path to run SimplusGT after installation.

% Author(s): Yitong Li, Yunjie Gu

function InstallSimplusGT()

clc

%% Generate msgbox
InstallMsg = 'Installing Simplus Grid Tool. This may take one minute...               ';
bar = waitbar(0.3,InstallMsg,'Name','Simplus');
bar.CloseRequestFcn = '';

%% Install
% Change the current folder
fprintf('Changing the current folder to SimplusGT folder...\n')
MfilePath = mfilename('fullpath');
[RootPath,~,~]  = fileparts(MfilePath);
cd(RootPath);

% Check the matlab version
fprintf('Checking the Matlab version...\n')
if hex2dec(version('-release')) < hex2dec('2016b')
    delete(bar);
    StopInstall('Error: Please use Matlab 2016a or later version!');
    return;
end

% Check if a previous version of SimplusGT has been installed
fprintf('Checking if the toolbox has been installed before...\n')
if exist('SimplusGT')~=0
    delete(bar);
    StopInstall(['Error: SimplusGT has been installed on this PC before or SimplusGT lib is opened. Please uninstall the old version first or close the SimplusGT lib.']);
    return;
end

% Add folder to path
fprintf('Installing...\n')
addpath(genpath([RootPath,'/Examples']));         	% Add "Examples" folder
addpath(genpath([RootPath,'/Library']));          	% Add "Library" folder
addpath(genpath([RootPath,'/Documentations']));   	% Add "Documentations" folder
addpath(genpath([RootPath,'/App']));                % Add "App" folder
addpath(genpath(RootPath));                         % Add root path
savepath;

% Waitbar
waitbar(0.6,bar,InstallMsg,'Name','Simplus');

% Convert the lib to the required version of matlab
fprintf('Converting SimplusGT library to the required Matlab version, please wait a minute...\n')
warning('off','all')    % Turn off all warnings
load_system('SimplusGT_2016a.slx');
save_system('SimplusGT_2016a.slx','Library/SimplusGT.slx');
close_system('SimplusGT.slx');
warning('on','all')     % Turn on all warnings

%%
fprintf('Congratulations!\n')
fprintf('Simplus Grid Tool successfully installed!\n')
fprintf('Please run "UserMain.m" later in the root path for using SimplusGT.\n')
fprintf('Please read manuals in the "Documentations" folder for more information.\n')

%% Remove msgbox
beep;
waitbar(1,bar,InstallMsg,'Name','Simplus');
pause(1);
delete(bar);

msgbox('Simplus Grid Tool successfully installed!','Simplus');

end

function StopInstall(msg)
    beep;
    uiwait(msgbox(msg,'Simplus','warn','modal'));
end