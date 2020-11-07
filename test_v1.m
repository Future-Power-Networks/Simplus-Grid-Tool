clear all
close all
clc

SCR = 2;
n = 10;

Xf = 0.03;
Rf = 0.01;

x = sym('x');

Xg = solve((1/n*x + Rf)^2 + (x + Xf)^2 - (1/SCR^2) == 0,x);
index = find(Xg>0);
Xg = vpa(Xg(index),4);