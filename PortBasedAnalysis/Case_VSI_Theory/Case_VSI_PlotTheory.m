% This function calculates the theoratical results of the
% single-VSI-infinite-bus system.

% Author(s): Yitong Li

%%
clear all;
clc;
close all;

enable_save = 1;
enable_LineEMT = 1;     % 1-proposed; 0-conventional.

s = sym('s');
F0 = 60;
W0 = F0*2*pi;
fn = 0;

omega_p = logspace(-2,4,10e3)*2*pi;
omega_pn = [-flip(omega_p),omega_p];

ColorRGB();

V_dc_r = 2.5;
C_dc = 1.25;

%%
% Load saved data
if enable_LineEMT
% With line EMT dynamcis
GsysDSS_SCR1d3  = load('GsysDSS_SCR1d3.mat').GsysDSS;
GsysDSS_SCR1d5  = load('GsysDSS_SCR1d5.mat').GsysDSS;
GsysDSS_SCR1d7  = load('GsysDSS_SCR1d7.mat').GsysDSS;
GsysDSS_SCR2    = load('GsysDSS_SCR2.mat').GsysDSS;
GsysDSS_SCR2d5  = load('GsysDSS_SCR2d5.mat').GsysDSS;
else
% Without line EMT dynamics
GsysDSS_SCR1d3  = load('st_GsysDSS_SCR1d3.mat').st_GsysDSS;
GsysDSS_SCR1d5  = load('st_GsysDSS_SCR1d5.mat').st_GsysDSS;
GsysDSS_SCR1d7  = load('st_GsysDSS_SCR1d7.mat').st_GsysDSS;
GsysDSS_SCR2    = load('st_GsysDSS_SCR2.mat').st_GsysDSS;
GsysDSS_SCR2d5  = load('st_GsysDSS_SCR2d5.mat').st_GsysDSS;
end

% Notes:
% These SCR values consider both the grid impedance and the output filter
% of inverter. If considering the grid impedance only, the updated SCR
% values should be:
% 1.3 -> 1.35
% 1.5 -> 1.57
% 1.7 -> 1.79
% 2   -> 2.13
% 2.5 -> 2.71
% See 'CalSCR.m' for more details.

% Get min realization
GminSS{5} = minreal(GsysDSS_SCR1d3);
GminSS{4} = minreal(GsysDSS_SCR1d5);
GminSS{3} = minreal(GsysDSS_SCR1d7);
GminSS{2} = minreal(GsysDSS_SCR2);
GminSS{1} = minreal(GsysDSS_SCR2d5);

% Calculate poles
pole_sys{5}     = pole(GsysDSS_SCR1d3)/2/pi; 
pole_sys{4}     = pole(GsysDSS_SCR1d5)/2/pi;
pole_sys{3}     = pole(GsysDSS_SCR1d7)/2/pi; 
pole_sys{2}     = pole(GsysDSS_SCR2)/2/pi; 
pole_sys{1}     = pole(GsysDSS_SCR2d5)/2/pi;

% Convert ss to sym for plotting
for i = 1:length(GminSS)
    G_dc_ss{i} = -GminSS{i}(7,6)*V_dc_r;
    G_dc{i} = -ss2sym(GminSS{i}(7,6))*V_dc_r;
end

% Calculate dc-current coefficient
for i = 1:length(G_dc)
    G_C = 1/(C_dc*s);
    K{i} = 1/G_dc{i} - 1/G_C;
    Gop{i} = K{i}*G_C;
    Gcl{i} = G_C/(1+Gop{i}); 
end

% Calculate the ss-form torque coefficient
for i = 1:length(G_dc)
    G_C_ss = 1/(C_dc*tf('s'));
    K_ss{i} = G_dc_ss{i}^(-1) - G_C_ss^(-1);
    Gop_ss{i} = K_ss{i}*G_C_ss;
end

% Calculate poles of K
for i = 1:length(K_ss)
    pole_K{i} = pole(K_ss{i});
end

%% Plot poles

if 1

fn = fn + 1;
figure(fn)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.35]);
for i = 1:length(pole_sys)
    scatter(real(pole_sys{i}),imag(pole_sys{i}),'x','LineWidth',1.5); hold on; grid on;
end
xlabel({'Real Part (Hz)'})
ylabel({'Imaginary Part (Hz)'})
X_H = 10;
X_L = -40;
XTick = [-40,-30,-20,-10,0,10];
set(gca,'XLim',[X_L,X_H]);
set(gca,'XTick',XTick);
Y_H = 150;
Y_L = -Y_H;
YTick = [-150,-100,-50,0,50,100,150];
set(gca,'YLim',[Y_L,Y_H]);
set(gca,'YTick',YTick);

if enable_save
    if enable_LineEMT
        print(gcf,'Case_VSI_PoleLocus.png','-dpng','-r600');
    else
        print(gcf,'Case_VSI_PoleLocus_NoLineEMT.png','-dpng','-r600');
    end
end

end

%% Plot K in complex plane

if 1
    
fn = fn+1;
figure(fn)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.25 0.35]);
if enable_LineEMT
    plot_K(K{1},81.04*1j*2*pi); grid on; hold on;
    plot_K(K{2},71.57*1j*2*pi); grid on; hold on;
    plot_K(K{3},65.27*1j*2*pi); grid on; hold on;
    plot_K(K{4},60.74*1j*2*pi); grid on; hold on;
    plot_K(K{5},55.86*1j*2*pi); grid on; hold on;
    % Notes:
    % The frequency shift is required because, for an inverter, the
    % frequency of the osillation mode is not close to 0 Hz anymore.
    X_L = -200;
    X_H = 300;
    set(gca,'XLim',[X_L,X_H]);
    set(gca,'XTick',[-200,0,200,300]);
    Y_L = -600;
    Y_H = 200;
    set(gca,'YLim',[Y_L,Y_H]);
    if enable_save
      print(gcf,'Case_VSI_Complex_K.png','-dpng','-r600');
    end
else
    plot_K(K{1},6.289*1j*2*pi); grid on; hold on;
    plot_K(K{2},6.022*1j*2*pi); grid on; hold on;
    plot_K(K{3},5.723*1j*2*pi); grid on; hold on;
    plot_K(K{4},5.286*1j*2*pi); grid on; hold on;
    plot_K(K{5},4.643*1j*2*pi); grid on; hold on;
  	X_L = 0;
    X_H = 200;
    set(gca,'XLim',[X_L,X_H]);
    set(gca,'XTick',[0,100,200]);
    Y_L = -300;
    Y_H = 0;
    set(gca,'YLim',[Y_L,Y_H]);
    set(gca,'YTick',[-300,-200,-100,0]);
	if enable_save
        print(gcf,'Case_VSI_Complex_K_NoLineEMT_zoomed.png','-dpng','-r600');
    end
   	X_L = -200;
    X_H = 300;
    set(gca,'XLim',[X_L,X_H]);
    set(gca,'XTick',[-200,0,200,300]);
    Y_L = -600;
    Y_H = 200;
    set(gca,'YLim',[Y_L,Y_H]);
    set(gca,'YTick',[-600,-400,-200,0,200]);
    if enable_save
        print(gcf,'Case_VSI_Complex_K_NoLineEMT.png','-dpng','-r600');
    end
end

end

%% Plot poles of K
if 0
fn = fn + 1;
figure(fn)
for i = [1,2,3,4]
    scatter(real(pole_K{i}),imag(pole_K{i}),'x','LineWidth',1.5); hold on; grid on;
end
title('pole:K')
end

%% Plot Bode diagram

if 0
fn = fn+1;
figure(fn)
legend_vec = {};
for i = 1:length(K)
    bodec(K{i},1j*omega_p,2*pi);
    legend_vec{i} = num2str(i);
end
legend(legend_vec)
title('K')
end

if 0
fn = fn+1;
figure(fn)
legend_vec = {};
for i = 1:length(Gop)
    bodec(Gop{i},1j*omega_p,2*pi);
    legend_vec{i} = num2str(i);
end
legend(legend_vec)
title('Gop:sym')
end

if 0
fn = fn+1;
figure(fn)
bode(Gop_ss{1}); grid on; hold on;
bode(Gop_ss{2}); grid on; hold on;
bode(Gop_ss{3}); grid on; hold on;
legend('1','2','3')
title('Gop:ss')
end
