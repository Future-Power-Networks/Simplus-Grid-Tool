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

Xline = 0.05;
Rline = Xline/5;
Lline = Xline/Wbase;

Xload1 = 0.2;
Rload1 = 1;
Lload1 = Xload1/Wbase;

Rload2 = 1;

w_pll = 1000*2*pi;
w_tau = 5000*2*pi;

