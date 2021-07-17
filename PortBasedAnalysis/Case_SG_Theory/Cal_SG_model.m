% This function gets the linearized state space model of a synchronous
% generator.

% Author(s): Yitong Li

clear all
close all
clc

% States
i_d = sym('i_d');
i_q = sym('i_q');
w = sym('w');

% Input
v_d = sym('v_d');
v_q = sym('v_q');
Tm = sym('Tm');

% Constant
psi_f = sym('psi_f');
R = sym('R');
L = sym('L');
J = sym('J');
D = sym('D');

%% State space equations in load convention

convention = 2;

% Input
u = [v_d; v_q; Tm];

% State
x = [i_d; i_q; w];

% Auxiliary equations
if convention == 1      % Motor convention
psi_d = L*i_d;
psi_q = L*i_q - psi_f;
Te = psi_f * i_d;
elseif convention == 2  % Generator convention
psi_d = -L*i_d;
psi_q = - L*i_q - psi_f;
Te = psi_f * i_d;
end

% state equations
if convention == 1          % Motor convention
di_d = (v_d - R*i_d + w*psi_q)/L;
di_q = (v_q - R*i_q - w*psi_d)/L;
dw = (Te - Tm -D*w)/J;
elseif convention == 2      % Generator convention
di_d = -(v_d + R*i_d + w*psi_q)/L;
di_q = -(v_q + R*i_q - w*psi_d)/L;
dw = (Tm - Te - D*w)/J;
end

f_xu = [di_d;
        di_q;
        dw];

% Output equations
g_xu = [i_d;
        i_q;
        w];
    
% Jacobian
A = jacobian(f_xu,x);
B = jacobian(f_xu,u);
C = jacobian(g_xu,x);
D = jacobian(g_xu,u);

%% Transfer function model
s = sym('s');
inv_sI_A = inv(s*eye(length(A)) - A);
G = C*inv_sI_A*B + D;
Gmin = simplify(G);

%% Output
A
B
C
D
G
