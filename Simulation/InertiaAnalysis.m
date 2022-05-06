% This script analyzes the inertia relationship between a SG and a IBR.
% The theory is based on the "Rethinking Grid-Fomring and Grid-Following
% ..." Paper.

% Author(s): Yitong Li

clear all
clc
close all

D = [7.5, 150, 150]
J = [300, 0.1, 1]

for i = 1:length(D)
    wc = D(i)/(2*J(i))/2/pi
end