%
% Version of UserMain using json input
%
% Author(s) Rob Oldaker
%

% Read me please:
% Please read the comments in this file carefully, and use this file to run
% toolbox.

%% Tips
%
% Please read manuals in the "Documentations" folder if you want to know
% more details about this tool.
%
% Please ensure that the toolbox is installed first, by running
% "InstallSimplusGT.m" once.
%
% The toolbox defaultly saves the results into Workspace, prints the key
% results in Command Window, and plots key figures.
%
% For changing default user data, please use "UserData.xlsx". More examples
% can be found in "Examples" folder.

%% Clear matlab
clear all;  % Clear Matlab workspace
clc;        % Clear Matlab command window
close all;  % Close all figures, etc

%% Set user data
% "UserData.xlsx" and "UserData.json" contain the data of an example 4-bus
% generator-inverter-composite power system. Please feel free to change
% them. 
% 
% ".xlsx" is the excel file, and ".json" is the corresponding json
% file. Users can easily convert an Excel file to a json file by calling
% this function saved in SimplusGT/Toolbox folder:
% SimplusGT.Toolbox.ConvertExcelFile2JsonFile();

% InputData = SimplusGT.JsonDecoder('UserData.json');

%%
% Other example power systems (in "Examples" folder):
%
% Pure ac power system examples:
% InputData = SimplusGT.JsonDecoder('SgInfiniteBus.json');              % Single synchronous generator and infinite bus
% InputData = SimplusGT.JsonDecoder('GflInverterInfiniteBus.json');   	% Single grid-following inverter and infinite bus
% InputData = SimplusGT.JsonDecoder('GfmInverterInfiniteBus.json');   	% Single grid-forming inverter and infinite bus
% InputData = SimplusGT.JsonDecoder('IEEE_14Bus.json');
% InputData = SimplusGT.JsonDecoder('IEEE_30Bus.json');
% InputData = SimplusGT.JsonDecoder('IEEE_57Bus.json');
% InputData = SimplusGT.JsonDecoder('NETS_NYPS_68Bus.json');
%
% Pure dc power system examples:
% InputData = SimplusGT.JsonDecoder('GfdBuckInfiniteBus.json');         % Single grid-feeding buck converter and infinite bus
%
% Hybrid ac-dc power system examples:
% InputData = SimplusGT.JsonDecoder('Hybrid_test_v1.json');             % A 4-bus hybrid ac-dc system
%
% For synchronisation analysis test:
% InputData = SimplusGT.JsonDecoder('Test_68Bus_IBR.json');
InputData = SimplusGT.JsonDecoder('Test_68Bus_IBR_17.json');
% InputData = SimplusGT.JsonDecoder('Test_68Bus_IBR_17_14.json');

%% Run toolbox
SimplusGT.Toolbox.MainStruct();

%% Results available to users (saved in Workspace)
% GsysDSS;          % Whole-system port model (descriptor state space
                    % form). Notes: The elements of state, input, and
                    % output vectors are printed in the command window.
                    %
                    % A quick introduction of DSS modeling method:
                    % https://uk.mathworks.com/help/simulink/slref/descriptorstatespace.html
                    
% GsysSS;           % Whole-system port model (state space form).
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
                    % flowing from each bus to the active apparatus connected.

% ListPowerFlow_;   % Power flow result for active apparatus only by combing 
                    % the PQ load into the nodal admittance matrix.
                    
% pole_sys;         % Whole-system poles, or equivalently eigenvalues.

% mymodel_v1;       % This is the simulink model generated automatically 
                    % based on the user data.