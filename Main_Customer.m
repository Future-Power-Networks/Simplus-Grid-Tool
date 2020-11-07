% This function is for customer use

%% Call "Main"
Main();

%% Available system
% These parameters are available for customer use
GsysObj;            % System object
GsysDSS;            % Descriptor state space model with original state preserved
GminSS;             % Minimum realization of descritptor state space model
ListPowerFlow;      % Power flow results
ListPowerFlow_;     % Power flow results (for active device only by combing the load into the nodal admittance matrix)
pole_sys;           % Poles of system

%% Customer plot