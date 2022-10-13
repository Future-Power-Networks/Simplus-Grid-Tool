% To users:
% Please use this file to run toolbox.
%% Clear matlab
clear all;  % Clear matlab workspace
clc;        % Clear matlab command window
close all;  % Close all figures, etc

%% Set user data
% 14 bus small-signal test
% UserData = 'MultiInverter_StrongGridInstability_Data.xlsx';

% 14 bus large-sigal (transient) test
% UserData = 'MultiInverter_Transient_GFM_Stable_Data.xlsx';
% UserData = 'MultiInverter_Transient_GFM_Unstable_Data.xlsx';
% UserData = 'MultiInverter_Transient_GFL_Data.xlsx';     

% UserData = 'MultiInverter_Transient_GFL_Load_Data.xlsx';     
% UserData = 'MultiInverter_Transient_GFM_Load_Data.xlsx';  
% Notes:
% The stability of GFM is adjusted by changing the droop gain in excel form.
% The stability of GFL is adjusted by changing ki_pll in "GridFollowingVSI.m".

UserData = 'SingleGflTest.xlsx';

%% Run toolbox
SimplexPS.Toolbox.Main();