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
syms Lf Cf

% Droop control
syms vqr wf Dw Dv Pr Qr Vd0 W0

% Voltage controller
syms kpv kiv

% Current controller
tau = 1/500;

%%
% System states
w = sym('w');
delta = sym('delta');
vdr = sym('vdr');
vdi = sym('vdi');
vqi = sym('vqi');
id = sym('id');
iq = sym('iq');
igd = sym('igd');
igq = sym('igq');
vd = sym('vd');
vq = sym('vq');

state = [w;delta;vdr;vdi;vqi;id;iq;igd;igq;vd;vq];

%%
% Droop control
% (Pr - P)*Dw/(1+s/wf) = W0 + w;
% (Qr - Q)*Dv/(1+s/wf) = Vd0 + vdr;
dw = ((Pr-p)*Dw+W0-w)*wf;
dvdr = ((Qr-q)*Dv+Vd0-vdr)*wf;

% Filter capacitor
dvd = id/Cf;
dvq = iq/Cf;
% dvd = (id-ild + Wg*Cf*vq)/Cf;
% dvq = (iq-ilq - Wg*Cf*vd)/Cf;

% Voltage controller
dvdi = vdr - vd;
dvqi = vqr - vq;

idr = kpv*dvdi + kiv*vdi;
iqr = kpv*dvqi + kiv*vqi;
% idr = kp*dvdi + ki*vdi - Cf*W0*vq;
% iqr = kp*dvqi + ki*vqi + Cf*W0*vd;

% Current controller
did = (idr - id)/tau;
diq = (iqr - iq)/tau;

% Angle difference between inverter and inf bus
% s*delta = w - wg;
ddelta = w - wg;

% Power calculation
p = vd*igd + vq*igq;
q = vq*igd - vd*igq;



% Frame transformation
% vD + j*vQ = (vd + j*vq)*e^{j*delta}
vD = vd*cos(delta) - vq*sin(delta);
vQ = vd*sin(delta) + vq*cos(delta);

% Line inductor
% vD - vgD = s*Lg*igd - wg*Lg*igq + Rg*igd
% vQ - vgQ = wg*Lg*igd + s*Lg*igq + Rg*igq;
digd = (vD - vgD + wg*Lg*igq - Rg*igd)/Lg;
digq = (vQ - vgQ - wg*Lg*igd - Rg*igq)/Lg;



%% Calculate the state matrix
f_xu = [dw;ddelta;dvdr;dvdi;dvqi;did;diq;digd;digq;dvd;dvq];

Amat = jacobian(f_xu,state);

%% Set numerical number
Cf = 0.02/Wbase;
Lf = 0.05/Wbase;

wf = 2*pi*10;

kpv = 0.00853;
kiv = 0.3062;

Dw = 0.05*Wbase;
Dv = 0.05;

vd = 1.0957;
vq = 0;
P = 1;
Q = 0.2279;
igd = P/vd;
igq = -Q/vd;
delta = 25.2846/180*pi;

Xg = 0.5;
Lg = Xg/Wbase;
Rg = Xg/5;

%% Replace symbolic by numerical number

Amat = subs(Amat,'kpv',kpv);
Amat = subs(Amat,'kiv',kiv);

Amat = subs(Amat,'Dw',Dw);
Amat = subs(Amat,'Dv',Dv);

Amat = subs(Amat,'vd',vd);
Amat = subs(Amat,'vq',vq);
Amat = subs(Amat,'delta',delta);
Amat = subs(Amat,'igd',igd);
Amat = subs(Amat,'igq',igq);

Amat = subs(Amat,'wf',wf);

Amat = subs(Amat,'Cf',Cf);
Amat = subs(Amat,'Lf',Lf);

Amat = subs(Amat,'Rg',Rg);
Amat = subs(Amat,'Lg',Lg);

Amat = subs(Amat,'wg',Wbase);
Amat = subs(Amat,'W0',Wbase);

%% Calculate pole
Amat = double(Amat)

EigVec = eig(Amat);
EigVecHz = EigVec/(2*pi);

PlotPoleMap(EigVecHz,9999);

