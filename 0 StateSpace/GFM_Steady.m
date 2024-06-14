% This function conducts the state space analysis of a grid-forming
% inverter. 
%
% Author(s): Yitong Li
% 
% The analysis is done mainly in steady-state frame analysis, i.e., DQ
% frame rather than dq frame. But the results in both frame-based analysis
% should be equivalent.
%
% By test, there is still an error between steady-frame-based and
% swing-frame-based state space analysis. In other words, the results of
% this function are slightly different from the results of
% "GfmModel_Swing_All.m". The reason is still unknown.

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

% Droop control
syms vm wf Dw Dv Pr Qr V0 W0

% Voltage controller
syms kpv kiv

% Current controller
syms kpi kii

% Cross-decoupling gain
Fcdv = 0;
Fcdi = 0;

%% System states
% Droop controller
syms w delta

% Voltage controller
syms vdi vqi

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

% Power calculation
p = vd*igd + vq*igq;
q = vq*igd - vd*igq;

% Droop control
dw = ((Pr-p)*Dw + W0 - w)*wf;
% dvm = ((Qr-q)*Dv + V0 - vm)*wf;

% Angle difference between inverter and inf bus
% s*delta = w - wg;
ddelta = w - wg;
dethata = w;

% Voltage controller
% dvdi = vm - vd;
dvdi = vm - vd;
dvqi = 0 - vq;
idr = kpv*dvdi + kiv*vdi - Fcdv*Cf*Wbase*vq;
iqr = kpv*dvqi + kiv*vqi + Fcdv*Cf*Wbase*vd;

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
state = [vdi; vqi; idi; iqi; iD; iQ; vD; vQ; igD; igQ; w; delta];
f_xu = [dvdi; dvqi; didi; diqi; diD; diQ; dvD; dvQ; digD; digQ; dw; ddelta];

Amat = jacobian(f_xu,state);

%% Set numerical number
% Just be careful that the equilibrium may be calculated approximately
Cf = 0.02/Wbase;
Lf = 0.05/Wbase;
Rf = 0.05/5;

wf = 2*pi*10;

wv = 250*2*pi;
kpv = Cf*wv;
kiv = Cf*wv^2/4*50;

wi = 1000*2*pi;
kpi = Lf*wi;
kii = Lf*(wi^2)/4;

Dw = 0.05*Wbase/Sbase;
Dv = 0;

vd = 1;
vq = 0;
vD = vd;
vQ = vq;

P = 0.5;
Q = 0.2;

igd = P/vd;
igq = -Q/vd;
igD = igd;
igQ = igq;

id = igd;
iq = igq;
iD = id;
iQ = iq;

delta = 5/180*pi;

Xg = 0.3;
Lg = Xg/Wbase;
Rg = Xg/5;

vgd = 1;
vgq = 0;
vgD = 1;
vgQ = 0;

ed = vd;
eq = vq;

idi = ed/kii;
iqi = eq/kii;

vdi = id/kiv;
vqi = iq/kiv;

vm = vd;

%% Replace symbolic by numerical number

Amat = subs(Amat,'kpv',kpv);
Amat = subs(Amat,'kiv',kiv);

Amat = subs(Amat,'kpi',kpi);
Amat = subs(Amat,'kii',kii);

Amat = subs(Amat,'vdi',vdi);
Amat = subs(Amat,'vqi',vqi);
Amat = subs(Amat,'idi',idi);
Amat = subs(Amat,'iqi',iqi);
Amat = subs(Amat,'vm',vm);

Amat = subs(Amat,'delta',delta);

Amat = subs(Amat,'Dw',Dw);
Amat = subs(Amat,'Dv',Dv);

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

Amat = subs(Amat,'wf',wf);

Amat = subs(Amat,'Cf',Cf);
Amat = subs(Amat,'Lf',Lf);
Amat = subs(Amat,'Rf',Rf);

Amat = subs(Amat,'Rg',Rg);
Amat = subs(Amat,'Lg',Lg);

Amat = subs(Amat,'w',Wbase);
Amat = subs(Amat,'wg',Wbase);

%% Calculate pole
Amat = double(Amat)

EigVec = eig(Amat);
EigVecHz = EigVec/(2*pi);

ZoomInAxis = [-20,10,-60,60];
PlotPoleMap(EigVecHz,ZoomInAxis,9999);

