% This script analyzes the inertia relationship between a SG and a IBR.
% The theory is based on the "Rethinking Grid-Fomring and Grid-Following
% ..." Paper.

% Author(s): Yitong Li

clear all
clc
close all

W0 = 2*pi*60;

% For PLL PI analysis
% The PLL controller is: G_PLL = kp_pll + ki_pll/s

f_pll = 10;
w_pll = f_pll*2*pi;

kp_pll = w_pll;
ki_pll = w_pll^2/4;

H = 1/ki_pll*W0/2;
D = H*kp_pll;

% For PLL LPF
% The PLL controller is: G_PLL = kp_pll / (1 + s/w_tau)

f_pll = 2;
w_pll = f_pll*2*pi;

kp_pll = w_pll;
ki_pll = 0;

f_tau = 100;
w_tau = f_tau*2*pi;

H = 1/(kp_pll*w_tau)*W0/2;
D = 1/(kp_pll*w_tau)*w_tau*W0; % 1/kp*W0