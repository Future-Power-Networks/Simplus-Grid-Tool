clear all
clc
close all

%%
enable_save = 0;

%%

s = tf('s');

F0 = 60;
W0 = F0*2*pi;

% L filter impedance
Lf = 0.05/W0;
Rf = 0.01;
Zf = s*Lf + Rf;

% C filter impedance
Cf = 0.0002/W0;
Yc = s*Cf;
Zc = 1/Yc;

% Line impedance
Ll = 0.3/W0;
Rl = 0.03;

Zl = s*Ll + Rl;

%% Calculate bandwidth for Voltage node
pole_w = pole(1/(Zl+Lf+Rf));
wb_v = H_CalcBandwidth(pole_w);
fb_v = wb_v/2/pi;

%%
% PI control impedance
w_i = 1000*2*pi;
kp_i = Lf*w_i;
ki_i = Lf*w_i^2/4;
Z_PIi = kp_i + ki_i/(s-1i*W0);

% Current node impedance
Zi = Z_PIi+Zf;
% Zi = Zi*Zc/(Zi+Zc);       % Ignore the parallel capacitance

% Channel gain
G = Zi / (Zi + Zl);
G = minreal(G);
pole_w = pole(G);
pole_f = pole_w/2/pi;
pole_pu = pole_f/F0;
figure(1)
scatter(real(pole_pu),imag(pole_pu),'x','LineWidth',1.5); hold on; grid on;
if enable_save
save('H_pole_pu','pole_pu');
end

%% Calculate bandwidth for current node
for k = 1:length(pole_w)
    wb_i_(k) = H_CalcBandwidth(pole_w(k));
end
wb_i = min(wb_i_);
fb_i = wb_i/2/pi;

%%
if 1
w_i = linspace(500,1000,5)*2*pi;

pole_f = {};
pole_pu = {};
for k = 1:length(w_i)
    kp_i = Lf*w_i(k);
    ki_i = Lf*w_i(k)^2/4;
    Z_PIi = kp_i + ki_i/(s-1i*W0);          
        % The channel bandwidth of current node depends on inner loop (current control) rather than outer loop (PLL).
    G = (Z_PIi+Zf) / (Zf + Z_PIi + Zl);
    G = minreal(G);
    pole_f{k} = pole(G)/2/pi;
    pole_pu{k} = pole_f{k}/F0;
end

if enable_save
save('H_pole_pu_sweep','pole_pu');
end

figure(2)
for k = 1:length(w_i)
    scatter(real(pole_pu{k}),imag(pole_pu{k}),'x','LineWidth',1.5); hold on; grid on;
end

% Notes:
% By sweeping w_i, we could find that the larger of w_i, the larger of w_b.

end
