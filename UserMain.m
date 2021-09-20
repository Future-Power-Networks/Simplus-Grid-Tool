% To users:
% Please use this file to run toolbox.

%% Tips
%
% Please ensure that the toolbox is installed first, by running
% "InstallSimplexPS.m" once.
%
% The toolbox defaultly saves the results into Workspace, prints the key
% results in Command Window, and plots key figures.
%
% For changing default user data, please use "UserData.xlsx". More examples
% can be found in "Examples" folder.

%% Clear matlab
clear all;  % Clear matlab workspace
clc;        % Clear matlab command window
close all;  % Close all figures, etc

%% Set user data
% Default
% "UserData.xlsx" defaultly contains the data of a 4-bus
% generator-inverter-composite power system. Please feel free to change it.
% UserData = 'UserData.xlsx';

% Examples used in the paper:
% UserData = 'H_SingleSgInfiniteBus.xlsx';
% UserData = 'H_SingleGflInfiniteBus.xlsx';
% UserData = 'Gamma_SingleSgInfiniteBus_ForSim';

% K analysis
% UserData = 'K_68Bus_SG_Load';
% UserData = 'K_68Bus_SG_IBR_Load';
% UserData = 'K_68Bus_SG_IBR';
% UserData = 'K_68Bus_SG_IBR_17';
% UserData = 'K_68Bus_SG_IBR_6';

% UserData = 'K_68Bus_SG_IBR_Test';
% UserData = 'K_68Bus_SG_IBR_Test_v2';
% UserData = 'K_68Bus_SG_IBR_Test_v3';
% Notes:       
% The system is stable when:
                                        % 1- increasing the sampling frequency
                                        % 2- adding the LPF for PLL, which seems to play a very important role in transient stability. 
                                        
% 4bus test
UserData = '4Bus_v1';


%% Run toolbox
SimplexPS.Toolbox.Main();

%% Results available to users (saved in Workspace)
% GsysDSS;          % Whole-system port model (descriptor state space
                    % form). Notes: The elements of state, input, and
                    % output vectors are printed in the command window.
                    %
                    % A quick introduction of DSS modeling method:
                    % https://uk.mathworks.com/help/simulink/slref/descriptorstatespace.html
                    
% GminSS;           % Whole-system port model (state space form).
                    % Notes: This model is the minimum realization of
                    % GsysDSS, which keeps the same input and output as
                    % GsysDSS, but reduces the order of state.

% YsysDSS;          % Whole-system admittance model (descriptor state space
                    % form). Notes: This model is derived from GsysDSS by
                    % selecting the voltage and current ports only and
                    % removing other input and output ports.
                    
% ListPowerFlow;    % Power flow
                    % Notes: The result is in the form of
                    % | bus | P | Q | V | angle | omega |
                    % P and Q are in load convention, i.e., the P and Q
                    % flowing from each bus to the active device connected.

% ListPowerFlow_;   % Power flow result for active device only by combing 
                    % the PQ load into the nodal admittance matrix.
                    
% pole_sys;         % Whole-system poles, or equivalently eigenvalues.

% mymodel_v1;       % This is the simulink model generated automatically 
                    % based on the user data.

%% User function
% Users can write their own functions here to further deal with the data
% mentioned above.