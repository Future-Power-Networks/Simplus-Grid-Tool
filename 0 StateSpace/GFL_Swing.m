% Author(s): Yitong Li

clear all
% close all
clc

%% Base value
BaseValue();

%% Set parameters
% Grid voltage
syms vgD vgQ wg

% Line impedance
syms Lg Rg

% LC Filter
syms Lf Cf Rf

% PLL controller
syms kp_pll ki_pll W0

% Current controller
syms kpi kii idr iqr

%% System states
% PLL controller
syms delta vqi

% Current controller
syms idi iqi

% Passive component
syms id iq vd vq igd igq

%%
% Inverse frame transformation
% igd + j*igq = (igD + j*igQ)*e^{-j*delta}
vgd = vgD*cos(delta) + vgQ*sin(delta);
vgq = -vgD*sin(delta) + vgQ*cos(delta);

% PLL control
% Equations:
% w = w0 + vq*(kp_pll + ki_pll/s)
dvqi = vq;
w = (kp_pll*vq + ki_pll*vqi) + W0;

% Angle difference between inverter and inf bus
% s*delta = w - wg;
ddelta = w - wg;
% dtheta = w;

% Current controller
didi = idr - id;
diqi = iqr - iq;
ed = kpi*didi + kii*idi;
eq = kpi*diqi + kii*iqi;

% Inverter-side inductor
did = (ed - vd + w*Lf*iq - Rf*id)/Lf;
diq = (eq - vq - w*Lf*id - Rf*iq)/Lf;

% Filter capacitor
dvd = (id-igd + w*Cf*vq)/Cf;
dvq = (iq-igq - w*Cf*vd)/Cf;

% Grid-side inductor
digd = (vd - vgd + w*Lg*igq - Rg*igd)/Lg;
digq = (vq - vgq - w*Lg*igd - Rg*igq)/Lg;

%% Calculate the state matrix
state = [idi; iqi; id; iq; vd; vq; igd; igq; vqi; delta];
f_xu = [didi; diqi; did; diq; dvd; dvq; digd; digq; dvqi; ddelta];

Amat = jacobian(f_xu,state);

%% Set numerical number
% Just be careful that the equilibrium may be calculated approximately
Cf = 0.02/Wbase;
Lf = 0.05/Wbase;
Rf = 0.01;

wi = 500*2*pi;
kpi = Lf*wi;
kii = Lf*(wi^2)/4;

wpll = 10*2*pi;
kp_pll = wpll;
ki_pll = wpll^2/4;

vd = 1.0817;
vq = 0;
vgD = 1;
vgQ = 0;

P = 0.5;
Q = 0.2;
igd = P/vd;
igq = -Q/vd;
igD = igd;
igQ = igq;
id = igd;
iq = igq;

delta = 3.6268/180*pi;

Xg = 0.3;
Lg = Xg/Wbase;
Rg = Xg/5;

vqi = 0;

%% Replace symbolic by numerical number

Amat = subs(Amat,'W0',Wbase);
Amat = subs(Amat,'vqi',vqi);

Amat = subs(Amat,'kp_pll',kp_pll);
Amat = subs(Amat,'ki_pll',ki_pll);

Amat = subs(Amat,'kpi',kpi);
Amat = subs(Amat,'kii',kii);

Amat = subs(Amat,'vd',vd);
Amat = subs(Amat,'vq',vq);
Amat = subs(Amat,'delta',delta);

Amat = subs(Amat,'igd',igd);
Amat = subs(Amat,'igq',igq);
% Amat = subs(Amat,'igD',igD);
% Amat = subs(Amat,'igQ',igQ);

Amat = subs(Amat,'vgD',vgD);
Amat = subs(Amat,'vgQ',vgQ);

Amat = subs(Amat,'id',id);
Amat = subs(Amat,'iq',iq);

Amat = subs(Amat,'Cf',Cf);
Amat = subs(Amat,'Lf',Lf);
Amat = subs(Amat,'Rf',Rf);

Amat = subs(Amat,'w',Wbase);
Amat = subs(Amat,'wg',Wbase);

Amat = subs(Amat,'Rg',Rg);
Amat = subs(Amat,'Lg',Lg);

%% Plot
Amat = double(Amat)

EigVec = eig(Amat);
EigVecHz = EigVec/(2*pi);
ZoomInAxis = [-20,10,-60,60];
PlotPoleMap(EigVecHz,ZoomInAxis,100);
