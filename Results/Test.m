clear all
clc
close all


V = 1+1i*0;

M = [1i*V,1i*V;
     -1i*conj(V),-1i*conj(V)];

[S_M1,S_M2,S_M3]= svd(M)
[E_M1,E_M2,E_M3] = eig(M)

[~,S_S_M1,~] = svd(S_M1)
[~,S_E_M1,~] = svd(E_M1)

M_ = transpose(M)*M

