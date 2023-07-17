% This model is the modified IEEE 14-BUS model in the paper.

%% Clear matlab
clear all;  % Clear Matlab workspace
clc;        % Clear Matlab command window
close all;  % Close all figures, etc

%% Step 1: choose the avaliable case study.
CaseStudy=12;

switch CaseStudy
    case 1; UserData = 'IEEE_14Bus_Cyprus_original.xlsx';  % IEEE-14 original model
    case 2; UserData = 'IEEE_14Bus_Cyprus_modified.xlsx';   % IEEE-14 Detuned model
    case 3; UserData = 'IEEE_14Bus_Cyprus_modified2.xlsx'; % IEEE-14 Tuned model
    case 4; UserData = 'NETS_NYPS_68_original'; % 68 bus original model.
    case 5; UserData = 'NETS_NYPS_68_modified_detuned'; % 68 bus detuned model.
    case 6; UserData = 'NETS_NYPS_68_modified_tuned'; % 68 bus tuned model.
    case 7; UserData = 'NETS_NYPS_68_modified_tuned_inter_area'; % 68 bus tuned to stabilise the 0.65Hz interarea mode.
    case 8; UserData = 'IEEE_14Bus_Cyprus_modified_SSS2.xlsx'; % for small-signal strength
    case 9; UserData = '4bus_case_a.xlsx'; % for large-signal strength: whole-system strength a
    case 10; UserData = '4bus_case_b.xlsx'; % for large-signal strength: whole-system strength  b
    case 11; UserData = 'NETS_NYPS_68_modified2.xlsx';
    case 12; UserData = 'TwoSG.xlsx';
    case 13; UserData = '4bus_case_c.xlsx'; % for large-signal strength: strength to connect
    case 14; UserData = 'IEEE_14Bus_Cyprus_modified_July.xlsx';
        % 4bus_case5: SG+SG+GFM+ GFL; case6: SG+SG+GFL+GFL
end

%UserData = 'SgInfiniteBus.xlsx'; % for large-signal strength 
UserData = 'SgInfiniteBus_FO2.xlsx';

%% Step 2: Run toolbox, get solutions of the system.
tic
SimplusGT.Toolbox.Main();
%toc
%ModalAnalysisAPP;

%% Step 3: Modal analysis: participation and eigenvalue sensitivity
% Step 3.1: run 'ModalInitialise.m'. 
% Each time you run this code, it will generate an excel file in 
% 'Simplus-Grid-Tool\Examples\ParticipationAnalysis\', from which you
% will need to select the interesting modes, apparatus, d-q axes...
% If you are using a different matlab version (not 2018b) as me, then run it and select
% manually. After the first time, if the system remains unchanged, don't run it
% otherwise you will loose your selections.
%
% SimplusGT.Modal.ModalInitialise();

% Step 3.2: run the analysis codes. Figures will pop out.
% Numerical resutls of the Sensitivity can be found in 'MdSensResult'.
% 
% SimplusGT.Modal.ModalAnalysis();
toc

%% Additional: Plot a map and show the propagation of the oscillation (for IEEE-14)
% To draw the map figure as shown in my paper. But you need to run
% ModalAnalysis first.
%run map_draw.m

%% Additional: Plot bode-plot from previous saved data. (for IEEE-14)
%run bode_draw.m


%% Additional: Pole-map plot (zoomed)
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

%% Additional: System solution results available to users (saved in Workspace)
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