% This function uninstalls the SimplexPowerSystem toolbox.

% To be continued.

%%
clear all
clc
close all

%% Remove folder from path
rmpath(genpath(pwd));
restoredefaultpath;
savepath;
clc;

fprintf('SimplexPowerSystem is uninstalled successfully! \n')
fprintf('Many thanks for choosing SimplexPowerSystem! \n')