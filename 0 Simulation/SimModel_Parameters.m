% This function generates the paramters at workspace for simulink model.

% Author(s): Yitong Li

%%
clear all
close all
clc

%%
t0 = 2;
dt = 0.2;

%% Fundamental parameters
Fs = 10;             % (kHz)
Fs = Fs*1e3;
Ts = 1/Fs;

%% Base values
Wbase = 2*pi*50;    % (rad/s)
Vbase = 1;
Sbase = 1;
Ibase = Sbase/Vbase;
Zbase = Vbase/Ibase;
Ybase = 1/Zbase;

%% AC filter parameters
Lf = 0.05;
Cf = 0.02;
Lc = 1e-9;

%% Grid-forming inverter
% Droop
m = 5/100;              % (pu)
w_droop = 15*2*pi;

% Current loop
w_i_GFM = 500*2*pi;

% Voltage loop
w_v_GFM = 250 *2*pi;
Scale_ki_v = 20;

% Rated line impedance
X_rated_GFM = 0.3;
R_rated_GFM = X_rated_GFM/5;

%% Grid-following inverter
% PLL
w_pll = 15 *2*pi;   % (rad/s)
w_tau = 1e3 *2*pi;  % (rad/s)

% Current loop
w_i_GFL = 250 *2*pi;    % (rad/s)

% Rated line impedance
X_rated_GFL = 0.4;
R_rated_GFL = X_rated_GFL/5;

