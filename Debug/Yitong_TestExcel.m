clear all
clc
close all

% [t1,t2] = xlsread('UserData','Apparatus')
% t1 = xlsread('UserData_VB_Yue','Apparatus')

A= [1,1,1;
    2,2,2;
    3,3,3];

t1 = A(1:2:end,:)