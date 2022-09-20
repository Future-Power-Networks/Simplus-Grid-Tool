clear all
clc
close all

syms L w s

T = [1,j;1,-j];
Apn = [(s+j*w)*L,0;
     0,(s-j*w)*L];
Adq = inv(T)*Apn*T;
Adq = simplify(Adq)