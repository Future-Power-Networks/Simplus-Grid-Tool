% This script analyzes the inertia relationship between a SG and a IBR.
% The theory is based on the "Revisiting Grid-Fomring and Grid-Following
% ..." Paper.

% Author(s): Yitong Li

clear all
clc
close all

D = [7.5, 150, 150]
J = [300, 0.1, 3]

for i = 1:length(D)
    fc = D(i)/(2*J(i))/2/pi         % in Hz
    tau = 1/(fc*2*pi)               % in s
    % Be careful that tau = 1/(2*pi*fc) rather than 1/fc.
end