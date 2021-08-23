% This script analyzes the inertia relationship between a SG and a IBR.
% The theory is based on the "Rethinking Grid-Fomring and Grid-Following
% ..." Paper.

% Author(s): Yitong Li

clear all
clc
close all

W0 = 2*pi*50;

f_pll = 10;
w_pll = f_pll*2*pi;

kp_pll = w_pll;
ki_pll = w_pll^2/4;

H = 1/ki_pll*W0/2
D = H*kp_pll