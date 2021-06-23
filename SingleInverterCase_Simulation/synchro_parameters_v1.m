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

% Droop
m = 5/100;

w_droop = 10*2*pi;      % This value will also influence the stability, but why? Because this value should not influence the steady state operating point.

% PLL
w_pll = 10*2*pi;
w_tau = 1e3*2*pi;
% w_pll = 5/100*Wbase
% w_tau = 10*2*pi;