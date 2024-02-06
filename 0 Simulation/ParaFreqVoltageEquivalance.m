clear all
clc
close all

Fs = 100*1e3;
Ts = 1/Fs;

Fbase = 50;
Wbase = Fbase*2*pi;
Vbase = 1;
Sbase = 1;
Ibase = Sbase/Vbase;
Zbase = Vbase/Ibase;

X = 0.98;
L = X/Wbase;
R = X/5;

w_pll = 1000*2*pi;
w_tau = 5000*2*pi;

