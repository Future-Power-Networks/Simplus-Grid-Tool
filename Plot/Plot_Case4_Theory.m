% This function plots the theoratical results for the IEEE 14 bus system.

% Author(s): Yitong Li

close all;
clear all;
clc

%%
% Common
fn = 0;

s = sym('s');
F0 = 60;
W0 = F0*2*pi;

omega_p = logspace(-1,3,1e3)*2*pi;
omega_pn = [-flip(omega_p),omega_p];

%% Load saved data
% Load system
GsysDSS_Z_1by5  = load('GsysDSS_Z_1by5.mat').GsysDSS;
% GsysDSS_Z_3by5 = load('GsysDSS_Z_3by5.mat').GsysDSS;
GsysDSS_Z_5by5 = load('GsysDSS_Z_5by5.mat').GsysDSS;

% Calculate pole
pole_1by5 = pole(GsysDSS_Z_1by5)/2/pi;
% pole_3by5 = pole(GsysDSS_Z_3by5)/2/pi;
pole_5by5 = pole(GsysDSS_Z_5by5)/2/pi;

%% Plot pole locus
    
fn = fn+1;
figure(fn);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.16 0.52]);
lw = 1.2;

subplot(2,1,1)
scatter(real(pole_5by5),imag(pole_5by5),'x','LineWidth',lw); hold on; grid on;
% scatter(real(pole_3by5),imag(pole_3by5),'x','LineWidth',lw); hold on; grid on;
scatter(real(pole_1by5),imag(pole_1by5),'x','LineWidth',lw); hold on; grid on;
xlabel({'Real (Hz)'},'interpreter','latex')
ylabel({'Imaginary (Hz)'},'interpreter','latex')
xlim([-1000,0]);
ylim([-2500,2500]);
yticks([-2500,-1250,0,1250,2500])

subplot(2,1,2)
scatter(real(pole_5by5),imag(pole_5by5),'x','LineWidth',lw); hold on; grid on;
% scatter(real(pole_3by5),imag(pole_3by5),'x','LineWidth',lw); hold on; grid on;
scatter(real(pole_1by5),imag(pole_1by5),'x','LineWidth',lw); hold on; grid on;
% xlabel({'Real (Hz)'},'interpreter','latex')
% ylabel({'Imaginary (Hz)'},'interpreter','latex')
xlim([-20,5]);
xticks([-20,-10,0,5]);
ylim([-30,30]);
yticks([-30,-15,0,15,30]);

if 1
    print(gcf,'Case4_GFM_Z.png','-dpng','-r600');
end

