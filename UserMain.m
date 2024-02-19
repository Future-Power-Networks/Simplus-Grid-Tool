% The case studies for Impedance Margin Ratio: a New Metric for Small-Signal System Strength
% Authors: Yue Zhu (yue.zhu18@imperial.ac.uk), etc.

% Notes: 
% * Add all folders and subfolders in the Matlab path before running.
% * My matlab version: R2022b
% * Certain Matlab add-ons may be required.

clear;
clc;
close all;

% choose the avaliable case studies: 1-6.
CaseStudy=1; 
switch CaseStudy
    case 1; UserData = '2IBR_4_bus_case1.xlsx'; % IMR 4-bus: same well tuned
    case 2; UserData = '2IBR_4_bus_case2.xlsx'; % IMR 4-bus: based on 24, higher PLL bandwidth to create control interactions
    case 3; UserData = '2IBR_4_bus_case3.xlsx'; % IMR 4-bus: based on 25, lower current control bandwidth for IBR-2 to 
    case 4; UserData = 'IEEE68_good.xlsx'; % IMR 68-bus: good tuned
    case 5; UserData = 'IEEE68_bad.xlsx'; % IMR 68-bus: bad tuned
    case 6; UserData = 'IEEE68_verybad.xlsx'; % IMR 68-bus: very bad tuned.
end
tic
SimplusGT.Toolbox.Main();
if CaseStudy>=3 % heatmap is avaliable for 68-bus case only.
    run CriticalIMR.m 
end
%toc
ModalAnalysisAPP; % Open the newly designed Modal analysis APP. This is only required once and results can be refreshed inside the APP.
toc