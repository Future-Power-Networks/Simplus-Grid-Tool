clear all
close all
clc

%% Base value
Fbase = 50;
Wbase = 2*pi*Fbase;
Vbase = 1;
Sbase = 1;
Ibase = Sbase/Vbase;
Zbase = Vbase/Ibase;
Ybase = 1/Zbase;

%% Set parameters
% Grid voltage
syms vgD vgQ wg

% Line impedance
syms Lg Rg

% LC Filter
syms Lf Cf Rf

% Droop control
syms vqr wf Dw Dv Pr Qr V0 W0

% Voltage controller
syms kpv kiv

% Current controller
syms kpi kii

%% System states
% Droop controller
syms w delta
syms vm

% Voltage controller
syms vdi vqi

% Current controller
syms idi iqi

% Passive component
syms id iq vd vq igD igQ

%%
% Inverse frame transformation
% igd + j*igq = (igD + j*igQ)*e^{-j*delta}
igd = igD*cos(delta) + igQ*sin(delta);
igq = -igD*sin(delta) + igQ*cos(delta);

% Power calculation
p = vd*igd + vq*igq;
q = vq*igd - vd*igq;
% p = vd*id + vq*iq;

% Droop control
% w = W0 + (Pr - P)*Dw/(1+s/wf);
% vm = V0 + (Qr - Q)*Dv/(1+s/wf);
dw = ((Pr-p)*Dw + W0 - w)*wf;
dvm = ((Qr-q)*Dv + V0 - vm)*wf;

% Angle difference between inverter and inf bus
% s*delta = w - wg;
ddelta = w - wg;

% Voltage controller
dvdi = vm - vd;
dvqi = 0 - vq;
idr = kpv*dvdi + kiv*vdi;
iqr = kpv*dvqi + kiv*vqi;
% idr = kp*dvdi + ki*vdi - Cf*W0*vq;
% iqr = kp*dvqi + ki*vqi + Cf*W0*vd;

% Current controller
didi = idr - id;
diqi = iqr - iq;
ed = kpi*didi + kii*idi;
eq = kpi*diqi + kii*iqi;

% Inverter-side inductor
% ed - vd = s*Lf*id - w*Lf*iq + Rf*id
% eq - vq = w*Lf*id + s*Lf*iq + Rf*iq
did = (ed - vd + w*Lf*iq - Rf*id)/Lf;
diq = (eq - vq - w*Lf*id - Rf*iq)/Lf;

% Filter capacitor
dvd = (id-igd + w*Cf*vq)/Cf;
dvq = (iq-igq - w*Cf*vd)/Cf;

% Frame transformation
% vD + j*vQ = (vd + j*vq) * e^{j*delta}
vD = vd*cos(delta) - vq*sin(delta);
vQ = vd*sin(delta) + vq*cos(delta);

% Grid-side inductor
% vD - vgD = s*Lg*igD - wg*Lg*igQ + Rg*igD
% vQ - vgQ = wg*Lg*igD + s*Lg*igQ + Rg*igQ
digD = (vD - vgD + wg*Lg*igQ - Rg*igD)/Lg;
digQ = (vQ - vgQ - wg*Lg*igD - Rg*igQ)/Lg;

% Is this modeling right? Or should I model Lf, Cf also in wg?

%% Calculate the state matrix
state = [w; delta; vm; vdi; vqi; idi; iqi; id; iq; vd; vq; igD; igQ];
f_xu = [dw; ddelta; dvm; dvdi; dvqi; didi; diqi; did; diq; dvd; dvq; digD; digQ];

Amat = jacobian(f_xu,state);

%% Set numerical number
Cf = 0.02/Wbase;
Lf = 0.05/Wbase;
Rf = 0.01;

wf = 2*pi*10;

wv = 250*2*pi;
kpv = Cf*wv;
kiv = Cf*wv^2/4*50;

wi = 1000*2*pi;
kpi = Lf*wi;
kii = Lf*(wi^2)/4;

Dw = 0.05*Wbase/Sbase;
Dv = 0;

vd = 1.1509;
vq = 0;
P = 0.5;
Q = 0.5;
igd = P/vd;
igq = -Q/vd;
igD = igd;
igQ = igq;
id = igd;
iq = igq;

delta = 5.9874/180*pi;

Xg = 0.3;
Lg = Xg/Wbase;
Rg = Xg/5;

%% Replace symbolic by numerical number

Amat = subs(Amat,'kpv',kpv);
Amat = subs(Amat,'kiv',kiv);

Amat = subs(Amat,'kpi',kpi);
Amat = subs(Amat,'kii',kii);

Amat = subs(Amat,'Dw',Dw);
Amat = subs(Amat,'Dv',Dv);

Amat = subs(Amat,'vd',vd);
Amat = subs(Amat,'vq',vq);
Amat = subs(Amat,'delta',delta);

% Amat = subs(Amat,'igd',igd);
% Amat = subs(Amat,'igq',igq);
Amat = subs(Amat,'igD',igD);
Amat = subs(Amat,'igQ',igQ);

Amat = subs(Amat,'id',id);
Amat = subs(Amat,'iq',iq);

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

PlotPoleMap(EigVecHz,9999);

