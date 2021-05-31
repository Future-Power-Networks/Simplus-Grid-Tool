%% Contact
% Yue Zhu
% yue.zhu18@imperial.ac.uk

%% Running instructions
% 1. You need to add this folder and all subfolders to the path before running this program. 
% In this m file, select the case you would like to analyze, then click run.
% 2. The simulink model will be opened after the computation finished. You can
% disable it to save time.
% 3. The Greybox model will be saved in the variables 'MDlayer1', 'MDlayer2'
% and 'MDlayer3' as structures, with field names.
% 4. If you are using Matlab version before 2019a, you have to switch the 
% simulink library before run time-domain simulation. Go to folder 'Library',
% delete 'SimplexPS.slx', and rename 'SimplexPS_2015a.slx' into 'SimplexPS.slx'. 
% 5. If you hope to see time-domain results directly, just go to
% '\StepResponsePlot\StepPlot.m' to plot my saved data. The data are saved
% from Simulink scopes directly.

%% Attention
% 1) the phrases 'tic', 'toc' are used here to estimate the solution (computation) 
% time for the system analysis. For IEEE-14 bus system, the solution time is about 15 s.
% For NETS-NYPS 68 bus system, the solution time is about 50 s.

% 2) Bus 1 and Bus 16 are swapt in this NETS_NYPS_68bus model!
% Because Bus1 has to be the Swing bus for the toolbox. But the results are swapt
% back in the paper.

% 3) FS mode is a specially customized mode to save some time-domain
% simulation time: At the first 5 seconds, the damping torchque coefficient
% D of each generator is deliberately set as a high value to stabilize the
% system quickly, then the value was changed back after t=5 s.
% The fast mode will NOT affect any of the frequency-domain results.

% 4) If you change any parameters, you have to run UserMain.m before
% running any simulink model; You also have to run ModalPrerun.m to creat a
% new Model excel form, as the modes will change. 

%% Clear matlab
clear all;  % Clear matlab workspace
clc;        % Clear matlab command window
close all;  % Close all figures, etc

%% Main selections 
GreyboxCaseSel = 1; % 1 : detuned 68 bus system. 
                    % 2 : tuned following the grey-box approach. 
                    % 3 : tuned against the grey-box approach. 
                    % 4 : IEEE-14 bus system in the first version of the paper.
SimulinkAutoOpen = 0; % 0: not opening simulink model after calculations.
                      %1: open simulink model after calculations.
First_time = 0; 
% 1 enable; 0 disable.
% * If you are using Matlab 2019a (same as me), then select 0.
% * If not, and this is your first time to run these codes, select 1.
% *If selecting 1, an excel file named as 'ModalConfig.....xlsx' will be
% generated autmatically in 'NetworkConfigData' folder. Then you have to select
% modes and devices manully in that file before doing the Greybox analysis.
% Sheet-1 is for the classic state participation factor.
% Sheet-2 is for the Grey-box three-layer analyze.
% Sheet-3 is for enbaling configurations.
% After selection, run('SimplexPS.Modal.ModalAnalysis.m'); 

%% codes.       
tic                    
if GreyboxCaseSel==1
    Name_Netlist = 'NETS_NYPS_68bus_detuned_gov_FS';
    SimplexPS.Toolbox.Main();
    toc
    if SimulinkAutoOpen==1
        open Case_68_withGOV_FSmode.slx;
    end
elseif GreyboxCaseSel==2
    Name_Netlist = 'NETS_NYPS_68bus_tuned_gov_FS';
    SimplexPS.Toolbox.Main();
    toc
    if SimulinkAutoOpen==1
        open Case_68_withGOV_FSmode.slx;
    end
elseif GreyboxCaseSel==3
    Name_Netlist = 'NETS_NYPS_68bus_tunedAgainst_gov_FS';
    SimplexPS.Toolbox.Main();
    toc
    if SimulinkAutoOpen==1
        open Case_68_withGOV_FSmode_Tuned_Against.slx;
    end
elseif GreyboxCaseSel==4
    Name_Netlist = 'IEEE_14Bus_detuned';
    SimplexPS.Toolbox.Main();
    toc
else
    error('Please select a valid case.')
end

if First_time == 1
    run('SimplexPS.Modal.ModalPreRun.m');
else
    run('SimplexPS.Modal.ModalAnalysis.m');
end
    