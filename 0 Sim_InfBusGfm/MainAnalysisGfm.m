clear all
clc
close all

%%
% s = tf('s');
s = sym('s','real');
w = sym('s','real');

%%
Fbase = 60;
Wbase = 2*pi*Fbase;

Xf = 0.05;
Lf = Xf/Wbase;
Rf = Xf/5;
Cf = 0.02/Wbase;
Xg = 0.2;
Lg = Xg/Wbase;
Rg = Xg/5;

w_i = 500;
w_v = 250;

kp_i = w_i*Lf;
ki_i = w_i^2*Lf/4;
kp_v = w_v*Cf;
ki_v = w_v^2*Cf/4*100;

Gdel = 1;

%%
Z_PI_i = (kp_i + ki_i/s)*Gdel;

Zin = Z_PI_i + (s+1i*w)*Lf + Rf;

Gi = Z_PI_i/Zin;
Y_PI_v = (kp_v + ki_v/s)*Gi;

Ypr = Y_PI_v + (s+1i*w)*Cf;

Ztot = 1/(1/Zin+Ypr) + (s+1i*w)*Lg + Rg;

