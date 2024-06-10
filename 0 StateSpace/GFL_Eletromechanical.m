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

% CL Filter
syms Cf Rd

% PLL controller
syms kp_pll ki_pll W0

% Current controller
syms wi  idr iqr 

%% System states
% PLL controller
syms delta vqi w

% Passive component
syms id iq igd igq vcd vcq

%%
% Inverse frame transformation
% igd + j*igq = (igD + j*igQ)*e^{-j*delta}
vgd = vgD*cos(delta) + vgQ*sin(delta);
vgq = -vgD*sin(delta) + vgQ*cos(delta);

% Angle difference between inverter and inf bus
% s*delta = w - wg;
ddelta = w - wg;
% dtheta = w;

% Current controller
did=(idr-id)*wi;
diq=(iqr-iq)*wi;

% Filter capacitor
icd = id - igd;
icq = iq - igq;
dvcd = (icd + w*Cf*vcq)/Cf;
dvcq = (icq - w*Cf*vcd)/Cf;
vd = vcd + icd*Rd;
vq = vcq + icq*Rd;

% Grid-side inductor
% vD - vgD = s*Lg*igD - wg*Lg*igQ + Rg*igD
% vQ - vgQ = wg*Lg*igD + s*Lg*igQ + Rg*igQ
digd = (vd - vgd + w*Lg*igq - Rg*igd)/Lg;
digq = (vq - vgq - w*Lg*igd - Rg*igq)/Lg;

% PLL control
% Equations:
% w = w0 + vq*(kp_pll + ki_pll/s)
dvqi = vq;
w_ = (kp_pll*vq + ki_pll*vqi) + W0;
wc = 2*pi*1000;
dw = (w_-w)*wc;

%% Calculate the state matrix
state = [ id; iq; vcd; vcq; igd; igq; vqi; delta; w];
f_xu = [ did; diq; dvcd; dvcq;  digd; digq; dvqi; ddelta; dw];

Amat = jacobian(f_xu,state);

%% Set numerical number
Cf = 0.02/Wbase;
Rd = 1e-2;

wi = 500*2*pi;

wpll = 10*2*pi;
kp_pll = wpll;
ki_pll = wpll^2/4;

vd = 1.0817;
vq = 0;
vgD = 1;
vgQ = 0;
vcd = 1;
vcq = 0;

Pr=0.5;
Qr=0.2;

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

Amat = subs(Amat,'wi',wi);

Amat = subs(Amat,'kp_pll',kp_pll);
Amat = subs(Amat,'ki_pll',ki_pll);

Amat = subs(Amat,'vqi',vqi);

Amat = subs(Amat,'vd',vd);
Amat = subs(Amat,'vq',vq);
Amat = subs(Amat,'delta',delta);
Amat = subs(Amat,'vcd',vcd);
Amat = subs(Amat,'vcq',vcq);

Amat = subs(Amat,'igd',igd);
Amat = subs(Amat,'igq',igq);
% Amat = subs(Amat,'igD',igD);
% Amat = subs(Amat,'igQ',igQ);

Amat = subs(Amat,'vgD',vgD);
Amat = subs(Amat,'vgQ',vgQ);

Amat = subs(Amat,'id',id);
Amat = subs(Amat,'iq',iq);

Amat = subs(Amat,'Cf',Cf);
Amat = subs(Amat,'Rd',Rd);

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

