% To users:
% Please use and ONLY use this file to run toolbox.

%% Tips

% Please ensure that the toolbox is installed， by running
% "InstallSimplexPS.m" the first time.

% The toolbox defaultly prints the results in Matlab command window, saves
% the results into Matlab workspace, and prints figures.

% For changing default data, please change "UserData.xlsx". More examples
% can be found in the "Examples" folder.

%% Clear matlab
clear all;  % Clear matlab workspace
clc;        % Clear matlab command window
close all;  % Close all figures, etc

%% Set the data
% Default
Name_Netlist = 'UserData.xlsx';
% "UserData.xlsx" defaultly contains the data of a 4-bus
% generator-inverter-composite power system. Please feel free to change it.

% Other example systems are:
% Name_Netlist = 'PortCoupling_SingleSGInfiniteBus.xlsx';  	% Single-generator-infinite-bus system
% Name_Netlist = 'PortCoupling_SingleVSIInfiniteBus.xlsx';	% Single-inverter-infinite-bus system
% Name_Netlist = 'PortCoupling_IEEE14Bus.xlsx';            	% 14 bus system
% Name_Netlist = 'Duality_SingleInverterInfiniteBus';
% Name_Netlist = 'TestGridForming';

%% Run toolbox
SimplexPS.Toolbox.Main();

%% Results available to users
% GsysDSS;          % Whole-system model (descriptor state space form)
                    % Note: The elements of state, input, and output
                    % vectors are printed in the command window.
                    
% GminSS;           % Whole-system model (state space form)
                    % Note: This model is the minimum realization of
                    % GsysDSS, which keeps the same input and output as
                    % GsysDSS, but reduces the order of state.
                    
% ListPowerFlow;    % Power flow result
                    % Note: The result is in the form of
                    % | bus | P | Q | V | angle | omega |

% ListPowerFlow_;   % Power flow result only for active device by combing 
                    % the load into the nodal admittance matrix.
                    
% pole_sys;         % Whole-system poles, or equivalently eigenvalues.

% mymodel_v1;       % This is the simulink model generated automatically 
                    % based on "CustomerData.xlsx".

%% User function and plot
% Users can write their own functions here to further deal with the data
% mentioned above.