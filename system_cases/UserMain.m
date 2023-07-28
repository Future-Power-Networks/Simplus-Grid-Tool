% This model is the modified IEEE 14-BUS model in the paper.

%% Clear matlab
clear all;  % Clear Matlab workspace
clc;        % Clear Matlab command window
close all;  % Close all figures, etc

%% case study for Disspating Energy Flow (DEF) method

DEF_case=2;
if DEF_case == 1
    UserData = 'DEF_2SG.xlsx'; % forced oscillation in 2 generation system
elseif DEF_case == 2
    UserData = 'DEF_IEEE14_case1.xlsx'; % a 8.8 Hz self-excited mode caused by delibarately detuned IBR current loop: DEF is effective
elseif DEF_case == 3
    UserData = 'DEF_IEEE14_case2.xlsx'; % a 1.5 Hz self-excited mode caused by high power output of IBR with low system strength: DEF fails
elseif DEF_case == 4
    UserData = 'DEF_2IBR_yue.xlsx'; % forced oscillation in 2 IBR system
elseif DEF_case == 5
    UserData = 'DEF_2GFL_Eugenie.xlsx';  % forced oscillation in 2 IBR system, including Lingling Fan's method
elseif DEF_case == 0 
    % not applied
end
% note-1: IBR curent satuation is set as 1.5/1.5 pu in case 1 and 4, and is set
% as 1.4/0.1 pu in case 2 and 3.
% note-2: for case 4, user could change IBR type from 10 to 11 for a constant PQ
% note-3: DEF method is included in the simulink model, and the calculation
% of ModeShape is in ModeShape_IEEE.slx model, which uses the recorded idq
% data from the original model. The data is not archieved in Github so
% users should un DEF_IEEE_casex.slx first which will record and save the
% data needed for mode-shape calcuation.


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