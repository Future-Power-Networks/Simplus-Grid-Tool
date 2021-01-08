clear all
clc
close all

a1 = sym('a1','real');
b1 = sym('b1','real');

a2 = sym('a2','real');
b2 = sym('b2','real');

t1 = a1 + 1i*b1;
t2 = a2 + 1i*b2;

t3 = imag(t1/t2)
t4 = imag(t1)/imag(t2)