% This function calculates the approximate SCR for the VSI infinite bus
% system.

% Author(s): Yitong Li

clear all
close all
clc

SCR = 1.7;
n = 10;

Xf = 0.03;
Rf = 0.01;

x = sym('x');

Xg = solve((1/n*x + Rf)^2 + (x + Xf)^2 - (1/SCR^2) == 0,x);
index = find(Xg>0);

Xg = vpa(Xg(index),4)
Rg = Xg/n

% X/R = 10
% P = 0.5
% Xg      stability     SCR(Zline+Zfilter)  SCR(Zline)
% 0.367   Y             2.5                 2.71
% 0.467   Y             2                   2.13
% 0.555   Y             1.7                 1.79
% 0.633   N             1.5                 1.57
% 0.735   N             1.3                 1.35
% 0.964   N             1                   1.03
