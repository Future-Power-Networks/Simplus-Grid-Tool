% This function calculates the impedance and swing equations for a
% frequency droop grid-forming inverter.

% Author(s): Yitong Li

clear all
clc

%%
% Laplace operator
s = tf('s');

% Symbolic calculation can also be used, but in this case, we should use
% solve function to calculate the roots.
% s = sym('s','real');
% root = solve(S_it__,s)/2/pi;

% Base frequency
Wbase = 2*pi*50;

%% Default parameters
% Equilibrium
I = 0 + 1i*0;
% V = 0.707 + 1i*(-0.707);
V = 1 + 1i*0;
W = Wbase;

% Droop gain
m = 5/100;

% Droop filter
w_f = 5*2*pi;

% Line impedance
% X/R ratio is 5 defaultly
Xline = 0.2;
Rline = Xline/5;

% LCL filter
Lf = 0.05/W;
Rf = 0.05/5;
Cf = 0.02/W;
Lc = 0/W;
Rc = 0;

% Control delay
Gdel = 1;

% Current loop
w_i = 1000*2*pi;

% Voltage loop
w_v = 250*2*pi;

%% Impedance analysis
% PI controllers
kp_i = w_i*Lf;
ki_i = w_i^2*Lf/4;

kp_v = w_v*Cf;
ki_v = w_v^2*Cf/4 * 20;
% kp_v = 1/(16*wi*Lf);
% ki_v = 1/(4*Lf);

% LPF
LPF = 1/(1+s/w_f);

% Droop controller
m = m*Wbase;
G_FD = m*LPF;

% Line impedance
Lline = Xline/W;
Zline = (s+1i*W)*Lline + Rline;
Zline_m = [Zline, 0;
         0,     conj(Zline)];

% Inner impedance
Z_PIi = (kp_i + ki_i/s)*Gdel;
Z_inner = Z_PIi + (s+1i*W)*Lf + Rf;
Z_inner_m = [Z_inner, 0;
             0,       conj(Z_inner)];
Gi_cl = Z_PIi/Z_inner;  % Current loop gain

% Parallel impedance
Y_PIv = (kp_v + ki_v/s)*Gi_cl;
Y_parallel = Y_PIv + (s+1i*W)*Cf;
Y_parallel_m = [Y_parallel, 0;
                0,          conj(Y_parallel)];
Gv_cl = 1/(1/Z_inner + Y_parallel)*Y_PIv;       % Voltage loop gain

% Outer impedance
Z_outer = (s+1i*W)*Lc + Rc + Zline;
Z_outer_m = [Z_outer, 0;
             0,       conj(Z_outer)];
         
% Whole-system impedance model
G_iw = 1/2*G_FD*[1 1];
I0 = [1i*I;
     -1i*conj(I)];
V0 = [1i*V;
      -1i*conj(V)];                    
Z_temp_m = inv(Y_parallel_m + inv(Z_inner_m));
Z_temp_m_prime = Z_temp_m + (V0 - Z_temp_m*I0)*inv(s + G_iw*I0)*G_iw;
Ztot_prime = Z_temp_m_prime + Z_outer_m;
Ytot_prime = inv(Ztot_prime);

% Calculate pole
pole_sys = pole(minreal(Ytot_prime));
pole_sys = pole_sys/2/pi;

ZoomInAxis = [-20,10,-60,60];
PlotPoleMap(pole_sys,ZoomInAxis,9999);

% %% Swing analysis
% % Swing equation: current-angle
% Sit = s*1/G_FD + 1i/2*(I - conj(I));
% Z_FD = 1/(2*Sit)*[1i*V,        1i*V;
%                   -1i*conj(V), -1i*conj(V)];
% % Sit_prime = Ztot_prime/Sit;           
% Sit_prime = Ztot_prime*inv(Zline_m);
% 
% % Calculate root
% root = pole(minreal(1/Sit));
% root_prime = pole(minreal(1/Sit_prime));
% 
% % Convert rad/s to Hz
% root = root/2/pi;
% root_prime = root_prime/2/pi;