clear all
close all
clc

%%
Bf = 0.02;
w0 = 50*2*pi;
Cf = Bf/w0;

tau = 1/500;
Tm = 50*1e-3;

kp = 0.00853;
ki = 0.3062;
% kp = kp*10;
% ki = ki*1000;

Dw = 0.05*w0;
Dv = 0.05;

Rl = sym('Rl');
Ll = sym('Ll');
Pr = sym('Pr');
Qr = sym('Qr');
vd0 = sym('vd0');

vqr = sym('vqr');
vgD = sym('vgD');
vgQ = sym('vgQ');

% vqr = 0;
% vgD = 1;
% vgQ = 0;

%%
w = sym('w');
delta = sym('delta');
vdr = sym('vdr');
vdi = sym('vdi');
vqi = sym('vqi');
id = sym('id');
iq = sym('iq');
ild = sym('ild');
ilq = sym('ilq');
vd = sym('vd');
vq = sym('vq');

state = [w;delta;vdr;vdi;vqi;id;iq;ild;ilq;vd;vq];

%%
% Voltage controller
dvdi = vdr - vd;
dvqi = vqr - vq;

idr = kp*dvdi + ki*vdi;
iqr = kp*dvqi + ki*vqi;

% Current controller
did = (idr - id)/tau;
diq = (iqr - iq)/tau;

% Angle difference
% s*delta = w - w0;
ddelta = w - w0;

% Filter capacitor
dvd = id/Cf;
dvq = iq/Cf;
% dvd = (id-ild + w0*Cf*vq)/Cf;
% dvq = (iq-ilq - w0*Cf*vd)/Cf;

% Frame transformation
% vD + j*vQ = (vd + j*vq)*e^{j*delta}
vD = vd*cos(delta) - vq*sin(delta);
vQ = vd*sin(delta) + vq*cos(delta);

% Line inductor
% vD - vgD = s*Ll*ild - w0*Ll*ilq + Rl*ild
% vQ - vgQ = w0*Ll*ild + s*Ll*ilq + Rl*ilq;
dild = (vD - vgD + w0*Ll*ilq - Rl*ild)/Ll;
dilq = (vQ - vgQ - w0*Ll*ild - Rl*ilq)/Ll;

p = vd*ild + vq*ilq;
q = vq*ild - vd*ilq;

% (Pr - P)*Dw/(1+Tm*s) = w0 + w;
% (Qr - Q)*Dv/(1+Tm*s) = vd0 + vdr;
dw = ((Pr-p)*Dw+w0-w)/Tm;
dvdr = ((Qr-q)*Dv+vd0-vdr)/Tm;

%%
f_xu = [dw;ddelta;dvdr;dvdi;dvqi;did;diq;dild;dilq;dvd;dvq];

Amat = jacobian(f_xu,state)

%%
vd = 1.0957;
vq = 0;
p = 1;
q = 0.2279;
ild = p/vd;
ilq = -q/vd;
delta = 25.2846/180*pi;

Amat = subs(Amat,'vd',vd);
Amat = subs(Amat,'vq',vq);
Amat = subs(Amat,'delta',delta);
Amat = subs(Amat,'ild',ild);
Amat = subs(Amat,'ilq',ilq);

Ll = 0.4903/w0;
Rl = 0.4903/5;
AmatNew = Amat;
AmatNew = subs(AmatNew,'Rl',Rl);
AmatNew = subs(AmatNew,'Ll',Ll);

AmatNew = double(AmatNew)

EigVec = eig(AmatNew);
EigVecHz = EigVec/(2*pi);

PlotPoleMap(EigVecHz,9999);

%%
function PlotPoleMap(EigVecHz,FigN)
  
    figure(FigN);
    
    subplot(1,2,1)
    scatter(real(EigVecHz),imag(EigVecHz),'x','LineWidth',1.5); hold on; grid on;
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
    title('Global pole map');
    
	subplot(1,2,2)
    scatter(real(EigVecHz),imag(EigVecHz),'x','LineWidth',1.5); hold on; grid on;
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
    title('Zoomed pole map');
    axis([-80,20,-150,150]);
end
