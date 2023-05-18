%% Clear matlab
clear all;  % Clear matlab workspace
clc;        % Clear matlab command window
% close all;  % Close all figures, etc

%% Set user data
% K analysis
% UserData = 'K_68Bus_IBR_Load';
% UserData = 'K_68Bus_IBR';

% SG System
% UserData = 'K_NETS_NYPS_68Bus';
% UserData = 'K_NETS_NYPS_68Bus_test';
% UserData = 'K_NETS_NYPS_68Bus_LowDamping';
% UserData = 'K_NETS_NYPS_68Bus_HighDamping';

% Notes:       
% The system is stable when:
                                        % 1- increasing the sampling frequency
                                        % 2- adding the LPF for PLL, which seems to play a very important role in transient stability. 
                                        
Enable_ParticipationFactorAnalysis = 0; % 1/0
ImagMax = 60;                           % Hz, the upper limit of the selected mode
ImagMin = 10;                           % Hz, the lower limit of the selected mode
                  
% Synchronization of GFL test
% UserData = 'Test_SingleGflInfBus';

% Test AVC
% UserData = 'Test_AVC';
UserData = 'K_68Bus_IBR_Load_AVC';

%% Run toolbox
SimplexPS.Toolbox.Main();