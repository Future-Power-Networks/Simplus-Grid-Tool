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
% Default
% "UserData.xlsx" defaultly contains the data of a 4-bus
% generator-inverter-composite power system. Please feel free to change it.
% UserData = 'UserData.xlsm';

% Other example power systems (in "Examples" folder):
%
% Pure ac power system examples:
% UserData = 'SgInfiniteBus.xlsx';              % Single synchronous generator and infinite bus
% UserData = 'GflInverterInfiniteBus.xlsx';   	% Single grid-following inverter and infinite bus
% UserData = 'GfmInverterInfiniteBus.xlsx';   	% Single grid-forming inverter and infinite bus
% UserData = 'IEEE_14Bus.xlsx';
% UserData = 'IEEE_30Bus.xlsx';
% UserData = 'IEEE_57Bus.xlsx';
% UserData = 'NETS_NYPS_68Bus.xlsx';
%
% Pure dc power system examples:
% UserData = 'GfdBuckInfiniteBus.xlsx';         % Single grid-feeding buck converter and infinite bus
%
% Hybrid ac-dc power system examples:
% UserData = 'Hybrid_test_v1.xlsx';             % A 4-bus hybrid ac-dc system
%
% For synchronisation test
% UserData = 'Test_68Bus_NETS_NYPS';      % Default NETS_NYPS system
% UserData = 'Test_68Bus_IBR_Load';       % IBRs with passvie loads
% UserData = 'Test_68Bus_IBR';            % IBRs with active loads
UserData = 'Test_68Bus_IBR_17';         % IBR at node 17 is repaced by a SG
% UserData = 'Test_68Bus_IBR_17_14';      % 17, 14
% UserData = 'Test_68Bus_IBR_17_14_7';    % 17, 14, 7
% UserData = 'Test_2Bus';
% UserData = 'Test_3Bus';

%% Run toolbox
SimplusGT.Toolbox.Main();

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
                    
%% User function
% Users can write their own functions here to further deal with the data
% mentioned above.