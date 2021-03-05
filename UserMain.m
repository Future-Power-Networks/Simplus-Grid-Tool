% To users:
% Please use and ONLY use this file to run toolbox.

%% Tips
%
% Please ensure that the toolbox is installed， by running
% "InstallSimplexPS.m" the first time.
%
% The toolbox defaultly prints the results in Command Window, saves the
% results into Workspace, and plots figures.
%
% For changing default user data, please change "UserData.xlsx". More
% examples can be found in the "Examples" folder.

%% Clear matlab
clear all;  % Clear matlab workspace
clc;        % Clear matlab command window
close all;  % Close all figures, etc

%% Set user data
% Default
UserData = 'UserData.xlsx';
% "UserData.xlsx" defaultly contains the data of a 4-bus
% generator-inverter-composite power system. Please feel free to change it.

% Other example systems (in "Examples" folder):
% Pure ac power system examples:
% UserData = 'SgInfiniteBus.slsx';              % Single synchronous generator and infinite bus
% UserData = 'GflInverterInfiniteBus.slsx';   	% Single grid-following inverter and infinite bus
% UserData = 'GfmInverterInfiniteBus.slsx';   	% Single grid-forming inverter and infinite bus
% UserData = 'IEEE_14Bus.slsx';
% UserData = 'IEEE_30Bus.xlsx';
% UserData = 'IEEE_57Bus.xlsx';
% UserData = 'NETS_NYPS_68Bus.xlsx';

% Pure dc power system examples:
% UserData = 'GfdBuckInfiniteBus.xlsx';          % Single grid-feeding buck converter and infinite bus

% Hybrid ac-dc power system examples:
UserData = 'Hybrid_test.xlsx';
UserData = 'Hybrid_test_v1.xlsx';

% For debug
% UserData = 'UserData_dc.xlsx';

%% Run toolbox
SimplexPS.Toolbox.Main();

%% Results available to users
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
                    % keeping the voltage and current ports only and
                    % removing other inputs and outputs.
                    
% ListPowerFlow;    % Power flow result
                    % Notes: The result is in the form of
                    % | bus | P | Q | V | angle | omega |

% ListPowerFlow_;   % Power flow result only for active device by combing 
                    % the load into the nodal admittance matrix.
                    
% pole_sys;         % Whole-system poles, or equivalently eigenvalues.

% mymodel_v1;       % This is the simulink model generated automatically 
                    % based on the user data.

%% User function
% Users can write their own functions here to further deal with the data
% mentioned above.