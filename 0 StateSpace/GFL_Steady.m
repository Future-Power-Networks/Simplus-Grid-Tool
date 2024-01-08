% This function conducts the state space analysis of grid-following
% inverter.
%
% Author(s): Yitong Li
%
% The analysis is done in steady-frame-based, i.e., DQ frame.
%
% The results of this function are still slightly different from the
% results of "GflModel_Swing.m". The reason is still unknown.

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

% PLL
syms w0 kp_pll ki_pll

% Current controller
syms kpi kii idr iqr

% Cross-decoupling gain
Fcdv = 0;
Fcdi = 0;

%% System states
% PLL
syms delta vqi

% Current controller
syms idi iqi

% Passive component
syms iD iQ vD vQ igD igQ

%%
% Inverse frame transformation
% idq = iDQ *e^{-j*delta}
id = iD*cos(delta) + iQ*sin(delta);
iq = -iD*sin(delta) + iQ*cos(delta);
% igdq = igDQ *e^{-j*delta}
igd = igD*cos(delta) + igQ*sin(delta);
igq = -igD*sin(delta) + igQ*cos(delta);
% vdq = vDQ *e^{-j*delta}
vd = vD*cos(delta) + vQ*sin(delta);
vq = -vD*sin(delta) + vQ*cos(delta);

% PLL
dvqi = vq;
w = (kp_pll*vq + ki_pll*vqi) + w0;
% w_ = (kp_pll + ki_pll*vqi) + w0;
% dw = (w_ - w)*wc;

% Angle difference between inverter and inf bus
% s*delta = w - wg;
ddelta = w - wg;

% Current controller
didi = idr - id;
diqi = iqr - iq;
ed = kpi*didi + kii*idi - Fcdi*Lf*Wbase*iq;
eq = kpi*diqi + kii*iqi + Fcdi*Lf*Wbase*iq;

% Frame transformation
% eDQ = edq * e^{j*delta}
eD = ed*cos(delta) - eq*sin(delta);
eQ = ed*sin(delta) + eq*cos(delta);

% Inverter-side inductor
diD = (eD - vD + wg*Lf*iQ - Rf*iD)/Lf;
diQ = (eQ - vQ - wg*Lf*iD - Rf*iQ)/Lf;

% Filter capacitor
dvD = (iD-igD + wg*Cf*vQ)/Cf;
dvQ = (iQ-igQ - wg*Cf*vD)/Cf;

% Grid-side inductor
digD = (vD - vgD + wg*Lg*igQ - Rg*igD)/Lg;
digQ = (vQ - vgQ - wg*Lg*igD - Rg*igQ)/Lg;

%% Calculate the state matrix
state = [idi; iqi; iD; iQ; vD; vQ; igD; igQ; vqi; delta];
f_xu = [didi; diqi; diD; diQ; dvD; dvQ; digD; digQ; dvqi; ddelta];

Amat = jacobian(f_xu,state);

%% Set numerical number
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
vD = vd;
vQ = vq;

P = 0.5;
Q = 0.2;

igd = P/vd;
igq = -Q/vd;
igD = igd;
igQ = igq;
idr = igD;
iqr = igQ;

id = igd;
iq = igq;
iD = id;
iQ = iq;

delta = 3.6268/180*pi;

Xg = 0.3;
Lg = Xg/Wbase;
Rg = Xg/5;

vgd = 1;
vgq = 0;
vgD = vgd;
vgQ = vgq;

ed = vd;
eq = vq;

idi = ed/kii;
iqi = eq/kii;

%% Replace symbolic by numerical number

Amat = subs(Amat,'kpi',kpi);
Amat = subs(Amat,'kii',kii);

Amat = subs(Amat,'idi',idi);
Amat = subs(Amat,'iqi',iqi);

Amat = subs(Amat,'kp_pll',kp_pll);
Amat = subs(Amat,'ki_pll',ki_pll);
Amat = subs(Amat,'delta',delta);

Amat = subs(Amat,'id',id);
Amat = subs(Amat,'iq',iq);
Amat = subs(Amat,'iD',iD);
Amat = subs(Amat,'iQ',iQ);

Amat = subs(Amat,'vd',vd);
Amat = subs(Amat,'vq',vq);
Amat = subs(Amat,'vD',vD);
Amat = subs(Amat,'vQ',vQ);

Amat = subs(Amat,'igd',igd);
Amat = subs(Amat,'igq',igq);
Amat = subs(Amat,'igD',igD);
Amat = subs(Amat,'igQ',igQ);
Amat = subs(Amat,'idr',idr);
Amat = subs(Amat,'iqr',iqr);

Amat = subs(Amat,'Cf',Cf);
Amat = subs(Amat,'Lf',Lf);
Amat = subs(Amat,'Rf',Rf);

Amat = subs(Amat,'Rg',Rg);
Amat = subs(Amat,'Lg',Lg);

Amat = subs(Amat,'w0',Wbase);
Amat = subs(Amat,'w',Wbase);
Amat = subs(Amat,'wg',Wbase);

%% Calculate pole
Amat = double(Amat)

EigVec = eig(Amat);
EigVecHz = EigVec/(2*pi);

ZoomInAxis = [-20,10,-60,60];
PlotPoleMap(EigVecHz,ZoomInAxis,100);

