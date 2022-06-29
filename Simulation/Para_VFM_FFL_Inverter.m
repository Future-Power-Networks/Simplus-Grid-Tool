% Author(s): Yitong Li

%%
clear all
close all
clc

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
Lc = 0.02; 

%% Grid-forming inverter
% DC-link control
vdc_ref = 2.5;
% Cdc = 1.25;     % 2*0.1*vdc^2;
Cdc = 1.25;
Rdc = 0.2;
w_vdc = 5*2*pi;

% Droop
m = 5/100;              % (pu)
w_droop = 5*2*pi;
w_pll = 15*2*pi;
w_tau = 200*2*pi;

% Current loop
w_i = 600*2*pi;

% Voltage loop
w_v = 300*2*pi;

% Line impedance
Xline = 0.95;
Rline = Xline/5;

%% Theoratical analysis
% m = 5/100*Wbase = 2.5Hz

% J = 1/(m*w_droop)
% D = 1/m

% ki = m*w_droop = 1/J
% kp = ki/m = w_droop = D/J

% PLL bandwidth = droop bandwidth

% If input of the rotor is torque rather than power
% J_ = J/wr
% D_ = D/wr

% ki = 1/Cdc*w/vdc^2
% kp = R*w/vdc^2

% Update Cdc and Rdc by referring to LPF
Cdc = Wbase/(vdc_ref^2)/(m*Wbase*w_droop)
Rdc = w_droop/Wbase*vdc_ref^2
