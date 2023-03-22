clear all
clc
close all

syms L1 L2
E = diag([L1,0,L2,0,0]);
A = [0,1,0,0,0;
     -1,0,0,0,1;
     0,0,0,1,0;
     0,0,-1,0,1;
     0,-1,0,-1,0];
B = [0;0;0;0;1];
C = [0,0,0,0,1];
D = 0;

[A_,B_,C_,D_,Bd,Dd]= CallDss2Ss(A,B,C,D,E);

A22 = [8,-99,-2;
       5,0,0;
       1,0,0];
   
R = rank(A22);

N = NullRight(A22);
t1 = A22*N

N_ = NullLeft(A22);
N_ = N_.'
t2 = N_*A22

% [t1,t2,t3] = svd(A22)