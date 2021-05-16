% This function plots the theoratical results for the IEEE 14 bus system.

% Author(s): Yitong Li

close all;
clear all;
clc

%%
% Common
ColorRGB();
Color{1} = color1;    
Color{2} = color2;
Color{3} = color3;
Color{4} = color4;
Color{5} = color5;

BodeSize = [0.1 0.1 0.3 0.6];

fn = 0;

enable_save = 1;

s = sym('s');
F0 = 60;
W0 = F0*2*pi;

omega_p = logspace(-1,3,1e3)*2*pi;
omega_pn = [-flip(omega_p),omega_p];

%% Load saved data
GsysDSS_5Hz  = load('GsysDSS_VSI5Hz.mat').GsysDSS;
GsysDSS_10Hz = load('GsysDSS_VSI10Hz.mat').GsysDSS;
GsysDSS_15Hz = load('GsysDSS_VSI15Hz.mat').GsysDSS;
GsysDSS_20Hz = load('GsysDSS_VSI20Hz.mat').GsysDSS;
GsysDSS_25Hz = load('GsysDSS_VSI25Hz.mat').GsysDSS;
GsysDSS_30Hz = load('GsysDSS_VSI30Hz.mat').GsysDSS;

%% Calculate poles
pole_5Hz = pole(GsysDSS_5Hz)/2/pi;
pole_10Hz = pole(GsysDSS_10Hz)/2/pi;
pole_15Hz = pole(GsysDSS_15Hz)/2/pi;
pole_20Hz = pole(GsysDSS_20Hz)/2/pi;
pole_25Hz = pole(GsysDSS_25Hz)/2/pi;
pole_30Hz = pole(GsysDSS_30Hz)/2/pi;

%% Get min realization
Gmin{1} = minreal(GsysDSS_5Hz);
Gmin{2} = minreal(GsysDSS_10Hz);
Gmin{3} = minreal(GsysDSS_15Hz);
Gmin{4} = minreal(GsysDSS_20Hz);
Gmin{5} = minreal(GsysDSS_25Hz);
Gmin{6} = minreal(GsysDSS_30Hz);

%% Get transfer functions
Ind_SG   = [1, 2, 6];
Ind_VSI  = [3, 8];
Port_w   = [3, 8, 13, 0, 0, 22, 0, 29];
Port_T_m = [3, 7, 0, 0, 0, 19];

J{1} = 3.5*2/W0^2;
J{2} = 0.35*2/W0^2;
J{3} = 3.5*2/W0^2;
D{1} = 1.5/W0^2;
D{2} = 1/W0^2;
D{3} = 1.5/W0^2;

for j = 1:length(Gmin)
    for i = 1:length(Ind_SG)

        n = Ind_SG(i);

        G_Tw{i,j} = -Gmin{j}(Port_w(n),Port_T_m(n));    % The SG model is in load convention.

        G_J{i} = 1/(J{i}*s);
        G_Tw{i,j} = ss2sym(G_Tw{i,j});
        K{i,j} = inv(G_Tw{i,j}) - inv(G_J{i});
        Gop_Tw{i,j} = K{i,j}*G_J{i};
        Gcl_Tw{i,j} = G_J{i}/(1+Gop_Tw{i,j});

    end
end

%% Plot pole locus

if 0
    
fn = fn+1;
figure(fn);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.8]);

subplot(2,1,1)
scatter(real(pole_5Hz),imag(pole_5Hz),'x','LineWidth',1.5); hold on; grid on;
scatter(real(pole_10Hz),imag(pole_10Hz),'x','LineWidth',1.5); hold on; grid on;
scatter(real(pole_15Hz),imag(pole_15Hz),'x','LineWidth',1.5); hold on; grid on;
scatter(real(pole_20Hz),imag(pole_20Hz),'x','LineWidth',1.5); hold on; grid on;
scatter(real(pole_25Hz),imag(pole_25Hz),'x','LineWidth',1.5); hold on; grid on;
xlabel({'Real Part (Hz)'})
ylabel({'Imaginary Part (Hz)'})
axis([-100,20,-200,200]);
set(gca,'XTick',[-100,-80,-60,-40,-20,0,20]);
set(gca,'YTick',[-200,-100,0,100,200]);

subplot(2,1,2)
scatter(real(pole_5Hz),imag(pole_5Hz),'x','LineWidth',1.5); hold on; grid on;
scatter(real(pole_10Hz),imag(pole_10Hz),'x','LineWidth',1.5); hold on; grid on;
scatter(real(pole_15Hz),imag(pole_15Hz),'x','LineWidth',1.5); hold on; grid on;
scatter(real(pole_20Hz),imag(pole_20Hz),'x','LineWidth',1.5); hold on; grid on;
scatter(real(pole_25Hz),imag(pole_25Hz),'x','LineWidth',1.5); hold on; grid on;
xlabel({'Real Part (Hz)'})
ylabel({'Imaginary Part (Hz)'})
axis([-0.15,0.05,-15,15]);
set(gca,'XTick',[-0.15,-0.1,-0.05,0,0.05]);
set(gca,'YTick',[-15,-10,-5,0,5,10,15]);
% set(gcf,'units','normalized','outerposition',[0.1 0.1 0.33 0.4]);

if enable_save
    print(gcf,'Case_Composite_PoleLocus.png','-dpng','-r600');
end

end

%% Transfer function from T to w
if 0

 	X_H = 2;
    X_L = 0;
    X_D = (X_H-X_L);
    XTick = logspace(X_L,X_H,X_D+1);
    XLim = [10^X_L,10^X_H];
    
 	Ym_H = 6;
    Ym_L = 0;
    Ym_D = (Ym_H - Ym_L);
    YTick_m = logspace(Ym_L,Ym_H,Ym_D+1);
    YLim_m = [10^Ym_L,10^Ym_H];

    Yp_H = 180;
    Yp_L = -180;
    Yp_D = (Yp_H - Yp_L)/90;
    YTick_p = linspace(Yp_L,Yp_H,Yp_D+1);
    YLim_p = [Yp_L,Yp_H];
    
    fn = fn+1;
    figure(fn);
    set(gcf,'units','normalized','outerposition',BodeSize);
   	for i=1:length(Ind_SG)
        bodec(G_Tw{i,5},1j*omega_p,2*pi,'PhaseOn',1,'Color',Color{i},'LineWidth',1.5);
    end
    legend('SG1','SG2','SG6','Location','Northeast')
  	f1 = subplot(2,1,1);
    f2 = subplot(2,1,2);
    set(f1,'XLim',XLim);
    set(f2,'XLim',XLim);
    set(f1,'XTick',XTick);
    set(f2,'XTick',XTick);
  	set(f1,'YLim',YLim_m);
    set(f1,'YTick',YTick_m);
  	set(f2,'YLim',YLim_p);
    set(f2,'YTick',YTick_p);
    
    % f1.XLabel.String = 'Frequency (Hz)';
    f1.YLabel.String = 'Magnitude';
    f2.XLabel.String = 'Positive Frequency (Hz)';
    f2.YLabel.String = 'Phase (Degree)';
    
    if enable_save
        print(gcf,'Case_Composite_Bode_Tw_25Hz.png','-dpng','-r600');
    end
    
end
    
%% Torque coefficient
if 1
    
	X_H = 2;
    X_L = 0;
    X_D = (X_H-X_L);
    XTick = logspace(X_L,X_H,X_D+1);
    XLim = [10^X_L,10^X_H];
    
    Ym_H = -1;
    Ym_L = -6;
    Ym_D = (Ym_H - Ym_L);
    YTick_m = logspace(Ym_L,Ym_H,Ym_D+1);
    YLim_m = [10^Ym_L,10^Ym_H];

    Yp_H = 180;
    Yp_L = -180;
    Yp_D = (Yp_H - Yp_L)/90;
    YTick_p = linspace(Yp_L,Yp_H,Yp_D+1);
    YLim_p = [Yp_L,Yp_H];
    
    % ### 5Hz
    fn = fn+1;
    figure(fn)
    set(gcf,'units','normalized','outerposition',BodeSize);
    for i=1:length(Ind_SG)
        bodec(K{i,1},1j*omega_p,2*pi,'PhaseOn',1,'Color',Color{i},'LineWidth',1.5);
    end
    % legend({'$K_1$','$K_2$','$K_6$'},'interpreter','latex','Location','Northeast')
    % legend({'$K_{1,N}$','$K_{2,N}$','$K_{6,N}$'},'interpreter','latex','Location','Northeast')
    %legend({'$K_{T1}$','$K_{T2}$','$K_{T6}$'},'interpreter','latex','Location','Northeast')
    legend({'$K_{T_{m1}}$','$K_{T_{m2}}$','$K_{T_{m6}}$'},'interpreter','latex','Location','Northeast')
   	f1 = subplot(2,1,1);
    f2 = subplot(2,1,2);
    
    % Magnified view
	set(f1,'XLim',[7,12]);
    set(f2,'XLim',[7,12]);
  	set(f1,'YLim',YLim_m);
    set(f1,'YTick',YTick_m);
  	set(f2,'YLim',[-180,0]);
    
  	f1.YLabel.String = 'Magnitude';
    f2.XLabel.String = 'Positive Frequency (Hz)';
    f2.YLabel.String = 'Phase (Degree)';
    
 	if enable_save
        print(gcf,'Case_Composite_Bode_K_5Hz_Zoomed.png','-dpng','-r600');
    end
    
    % Global view
    set(f1,'XLim',XLim);
    set(f2,'XLim',XLim);
    set(f1,'XTick',XTick);
    set(f2,'XTick',XTick);
  	set(f1,'YLim',YLim_m);
    set(f1,'YTick',YTick_m);
  	set(f2,'YLim',YLim_p);
    set(f2,'YTick',YTick_p);
    
  	f1.YLabel.String = 'Magnitude';
    f2.XLabel.String = 'Positive Frequency (Hz)';
    f2.YLabel.String = 'Phase (Degree)';
    
 	if enable_save
        print(gcf,'Case_Composite_Bode_K_5Hz.png','-dpng','-r600');
    end
    
    % ### 25 Hz
   	fn = fn+1;
    figure(fn)
    set(gcf,'units','normalized','outerposition',BodeSize);
    for i=1:length(Ind_SG)
        bodec(K{i,5},1j*omega_p,2*pi,'PhaseOn',1,'Color',Color{i},'LineWidth',1.5);
    end
    % legend({'$K_1$','$K_2$','$K_6$'},'interpreter','latex','Location','Northeast')
    % legend({'$K_{1,N}$','$K_{2,N}$','$K_{6,N}$'},'interpreter','latex','Location','Northeast')
    % legend({'$K_{T1}$','$K_{T2}$','$K_{T6}$'},'interpreter','latex','Location','Northeast')
    legend({'$K_{T_{m1}}$','$K_{T_{m2}}$','$K_{T_{m6}}$'},'interpreter','latex','Location','Northeast')
   	f1 = subplot(2,1,1);
    f2 = subplot(2,1,2);
    set(f1,'XLim',XLim);
    set(f2,'XLim',XLim);
    set(f1,'XTick',XTick);
    set(f2,'XTick',XTick);
 	set(f1,'YLim',YLim_m);
    set(f1,'YTick',YTick_m);
  	set(f2,'YLim',YLim_p);
    set(f2,'YTick',YTick_p);
    
  	f1.YLabel.String = 'Magnitude';
    f2.XLabel.String = 'Positive Frequency (Hz)';
    f2.YLabel.String = 'Phase (Degree)';
    
   	if enable_save
        print(gcf,'Case_Composite_Bode_K_25Hz.png','-dpng','-r600');
    end
    
end

%% Open loop gain
if 0 
    
    X_H = 2;
    X_L = 0;
    X_D = (X_H-X_L);
    XTick = logspace(X_L,X_H,X_D+1);
    XLim = [10^X_L,10^X_H];
    BodeSize = [0.1 0.1 0.35 0.5];
    
    Ym_H = 3;
    Ym_L = -3;
    Ym_D = (Ym_H - Ym_L);
    YTick_m = logspace(Ym_L,Ym_H,Ym_D+1);
    YLim_m = [10^Ym_L,10^Ym_H];

    Yp_H = 270;
    Yp_L = -270;
    Yp_D = (Yp_H - Yp_L)/90;
    YTick_p = linspace(Yp_L,Yp_H,Yp_D+1);
    YLim_p = [Yp_L,Yp_H];
    
    fn = fn+1;
    figure(fn)
    set(gcf,'units','normalized','outerposition',BodeSize);
    for i=1:length(Ind_SG)
        bodec(Gop_Tw{i,1},1j*omega_p,2*pi,'PhaseOn',1,'Color',Color{i},'LineWidth',1.5);
    end
    legend('SG1','SG2','SG6','Location','Northeast')
    f1 = subplot(2,1,1);
    f2 = subplot(2,1,2);
    set(f1,'XLim',XLim);
    set(f1,'XTick',XTick);
    set(f2,'XLim',XLim);
    set(f2,'XTick',XTick);
	set(f1,'YLim',YLim_m);
    set(f1,'YTick',YTick_m);
  	set(f2,'YLim',YLim_p);
    set(f2,'YTick',YTick_p);
    
    if enable_save
        print(gcf,'Case_Composite_Bode_Gop_5Hz.png','-dpng','-r600');
    end
    
   	fn = fn+1;
    figure(fn)
    set(gcf,'units','normalized','outerposition',BodeSize);
    for i=1:length(Ind_SG)
        bodec(Gop_Tw{i,5},1j*omega_p,2*pi,'PhaseOn',1,'Color',Color{i},'LineWidth',1.5);
    end
    legend('SG1','SG2','SG6','Location','Northeast')
   	f1 = subplot(2,1,1);
    f2 = subplot(2,1,2);
    set(f1,'XLim',XLim);
    set(f2,'XLim',XLim);
    set(f1,'XTick',XTick);
    set(f2,'XTick',XTick);
  	set(f2,'XTick',XTick);
	set(f1,'YLim',YLim_m);
    set(f1,'YTick',YTick_m);
  	set(f2,'YLim',YLim_p);
    set(f2,'YTick',YTick_p);
    
  	if enable_save
        print(gcf,'Case_Composite_Bode_Gop_25Hz.png','-dpng','-r600');
    end
    
end
