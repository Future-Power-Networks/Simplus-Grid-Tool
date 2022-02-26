% This model is the modified IEEE 14-BUS model in the paper.

%% Clear matlab
clear all;  % Clear Matlab workspace
clc;        % Clear Matlab command window
close all;  % Close all figures, etc

%% Set user data, select one excel file to form the system.
%UserData = 'IEEE_14Bus_Cyprus_original.xlsx';  % original model
UserData = 'IEEE_14Bus_Cyprus_modified.xlsx';   % Detuned model
% UserData = 'IEEE_14Bus_Cyprus_modified2.xlsx'; % Tuned model

%% Run toolbox, get solutions of the system.
SimplusGT.Toolbox.Main();

%% Modal analysis: participation and eigenvalue sensitivity
% Step-1: run 'ModalInitialise.m'. 
% Each time you run this code, it will generate an excel file in 
% 'Simplus-Grid-Tool\Examples\ParticipationAnalysis\', from which you can
% will need to select the interesting modes, apparatus, d-q axes...
% So if the system remains unchanged, don't run it.
%
% SimplusGT.Modal.ModalInitialise();

% Step-2: run the analysis codes. Figures will pop out.
% Numerical resutls of the Sensitivity can be found in 'MdSensResult'.
% 
%SimplusGT.Modal.ModalAnalysis();

%% Plot a map and show the propagation of the oscillation
%run map_draw.m

%% Plot bode-plot from previous saved data.
%run bode_draw.m


%% Pole-map plot (zoomed)
% this is used to plot the same zoomed pole-map as the original model of Cyprus
% in DigSLIENT Power Factory. The poles are nearly the same, proving that the 
% model is successfully transplanted to this tool box.

% figure(1);
% clf
% scatter(real(pole_sys)*2*pi,imag(pole_sys),'x','LineWidth',1.5); hold on; grid on;
% xlabel('Real Part (Hz)');
% ylabel('Imaginary Part (Hz)');
% title('Zoomed pole map');
% axis([-30,10,-3.5,3.5]);
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
                    % flowing from each bus to the active apparatus connected.

% ListPowerFlow_;   % Power flow result for active apparatus only by combing 
                    % the PQ load into the nodal admittance matrix.
                    
% pole_sys;         % Whole-system poles, or equivalently eigenvalues.

% mymodel_v1;       % This is the simulink model generated automatically 
                    % based on the user data.
                    
%% User function
% Users can write their own functions here to further deal with the data
% mentioned above.