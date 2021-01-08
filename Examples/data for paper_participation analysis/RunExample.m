% To customer:
% Please use and only use this file to run toolbox.

%% Tips

% Please ensure that the toolbox is installed (by using "InstallSimplexPS.m") the
% first time.

% For changing default data, please use "CustomerData.xlsx". (More examples
% can be found in the "Examples" folder)

% The toolbox defaultly prints the results in Matlab command window, saves
% the results into Matlab workspace, and prints figures.

%% Clear matlab
clear all;  % Clear matlab workspace
clc;        % Clear matlab command window
close all;  % Close all figures, etc

%% Run toolbox
Name_Netlist = '14Bus.xlsx';
% Name_Netlist = 'SingleSGInfiniteBus.xlsx';
% Name_Netlist = 'SingleVSIInfiniteBus.xlsx';
% Name_Netlist = '14Bus.xlsx';

SimplexPS.Toolbox.Main(); 	% This function runs toolbox.

%% Custormer available results
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

%% Customer function and plot
% Custormer can write their functions here to further deal with the data
% mentioned above.