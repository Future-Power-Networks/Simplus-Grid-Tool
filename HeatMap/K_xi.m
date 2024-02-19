% This function analyzes the xi

clear all
clc
close all

mfile_name = mfilename('fullpath');
[RootPath,~,~]  = fileparts(mfile_name);
cd(RootPath);

%% Enables
Enable_SaveFigure = 1;

%% Load data
xi{1}   = load('K_SG_Load_xi_diag').xi_diag;
xi{2}   = load('K_SG_IBR_Load_xi_diag').xi_diag;
xi{3}   = load('K_SG_IBR_xi_diag').xi_diag;
xi{4}   = load('K_SG_IBR_17_xi_diag').xi_diag;

%% Plot
for i = 1:4
    [~,index] = min(real(xi{i}));
    xi_min(i) = xi{i}(index);
end

plot(real(xi_min)); grid on; hold on;
