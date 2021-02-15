% This file uninstalls the SimplexPowerSystem toolbox.

% Author(s): Yitong Li

%%
clear all
clc
close all

%% 
% Remove folder from path
rmpath(genpath(pwd));
restoredefaultpath;
savepath;
clc;

fprintf('SimplexPowerSystem is uninstalled! \n')
fprintf('Many thanks for using SimplexPowerSystem! \n')
fprintf('Please delete the files manually for fully removing this toolbox from your PC. \n')