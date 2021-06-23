% This function generates the paramters at workspace for simulink model.

% Author(s): Yitong Li

clear all
close all
clc

Fs = 10;             % (kHz)
Fs = Fs*1e3;
Ts = 1/Fs;

Wbase = 2*pi*50;    % (rad/s)
Vbase = 1;
Sbase = 1;
Ibase = Sbase/Vbase;
Zbase = Vbase/Ibase;
Ybase = 1/Zbase;

%% Parameters
% PLL
w_pll = 15*2*pi;
w_tau = 1e3*2*pi;

w_i = 250*2*pi;     % (rad/s)

Lf = 0.05;
Cf = 0.02;
Lc = 1e-9;