% This function deals with the theoratical results of the
% single-SG-infinite-bus system.

% Author(s): Yitong Li

%%
clear all;
close all;

enable_save = 1;

s = sym('s');
F0 = 60;
W0 = F0*2*pi;
fn = 0;

omega_p = logspace(-2,4,10e3)*2*pi;
omega_pn = [-flip(omega_p),omega_p];

ColorRGB();

%%
% Load saved data
% GsysDSS_J3d5_D0d1_R0    = load('GsysDSS_J3d5_D0d1_R0.mat').GsysDSS_J3d5_D0d1_R0;
% GsysDSS_J3d5_D0d1_R0d05 = load('GsysDSS_J3d5_D0d1_R0d05.mat').GsysDSS_J3d5_D0d1_R0d05;
% GsysDSS_J3d5_D0d1_R0d1  = load('GsysDSS_J3d5_D0d1_R0d1.mat').GsysDSS_J3d5_D0d1_R0d1;
% GsysDSS_J3d5_D0d1_R0d2  = load('GsysDSS_J3d5_D0d1_R0d2.mat').GsysDSS_J3d5_D0d1_R0d2;
% GsysDSS_J3d5_D0d1_R0d3  = load('GsysDSS_J3d5_D0d1_R0d3.mat').GsysDSS_J3d5_D0d1_R0d3;

% GsysDSS_J3d5_D0d1_R0    = load('GsysDSS_psif1_D0d1_R0.mat').GsysDSS;
% GsysDSS_J3d5_D0d1_R0d05 = load('GsysDSS_psif1_D0d1_R0d05.mat').GsysDSS;
% GsysDSS_J3d5_D0d1_R0d1  = load('GsysDSS_psif1_D0d1_R0d1.mat').GsysDSS;
% GsysDSS_J3d5_D0d1_R0d2  = load('GsysDSS_psif1_D0d1_R0d2.mat').GsysDSS;
% GsysDSS_J3d5_D0d1_R0d3  = load('GsysDSS_psif1_D0d1_R0d3.mat').GsysDSS;

GsysDSS_J3d5_D0d1_R0    = load('GsysDSS_psif1_D0d2_R0.mat').GsysDSS;
GsysDSS_J3d5_D0d1_R0d05 = load('GsysDSS_psif1_D0d2_R0d05.mat').GsysDSS;
GsysDSS_J3d5_D0d1_R0d1  = load('GsysDSS_psif1_D0d2_R0d1.mat').GsysDSS;
GsysDSS_J3d5_D0d1_R0d2  = load('GsysDSS_psif1_D0d2_R0d2.mat').GsysDSS;
GsysDSS_J3d5_D0d1_R0d3  = load('GsysDSS_psif1_D0d2_R0d3.mat').GsysDSS;

% Get min realization
GminSS{1} = minreal(GsysDSS_J3d5_D0d1_R0);
GminSS{2} = minreal(GsysDSS_J3d5_D0d1_R0d05);
GminSS{3} = minreal(GsysDSS_J3d5_D0d1_R0d1);
GminSS{4} = minreal(GsysDSS_J3d5_D0d1_R0d2);
GminSS{5} = minreal(GsysDSS_J3d5_D0d1_R0d3);

% Calculate poles
pole_sys{1}     = pole(GsysDSS_J3d5_D0d1_R0)/2/pi; 
pole_sys{2}     = pole(GsysDSS_J3d5_D0d1_R0d05)/2/pi; 
pole_sys{3}     = pole(GsysDSS_J3d5_D0d1_R0d1)/2/pi;
pole_sys{4}     = pole(GsysDSS_J3d5_D0d1_R0d2)/2/pi; 
pole_sys{5}     = pole(GsysDSS_J3d5_D0d1_R0d3)/2/pi; 

% Convert ss to sym for plotting
for i = 1:length(GminSS)
    % Get the element
    G_Tw{i} = -ss2sym(GminSS{i}(6,5));
    G_vi{i} = -ss2sym(GminSS{i}(4:5,3:4));
    G_vw{i} = -ss2sym(GminSS{i}(6,3:4));
    G_Ti{i} = -ss2sym(GminSS{i}(4:5,5));
    
    % Convert matrix to complex form
    Tj = [1,1j;1,-1j];
    G_vi_c{i} = Tj*G_vi{i}*Tj^(-1);
    G_vw_c{i} = G_vw{i}*Tj^(-1);
    G_Ti_c{i} = Tj*G_Ti{i};
    
    % Get the ss form
    G_Tw_ss{i} = -GminSS{i}(6,5);
end

% Calculate torque coefficient
for i = 1:length(G_Tw)
    J = 3.5*2/W0^2;

    G = 1/(J*s);
    K{i} = 1/G_Tw{i} - 1/G;
    Gop{i} = K{i}/J/s;
    Gcl{i} = G/(1+Gop{i}); 
end

% Calculate the ss-form torque coefficient
for i = 1:length(G_Tw)
    J = 3.5*2/W0^2;
    G_ss = 1/(J*tf('s'));
    
    K_ss{i} = G_Tw_ss{i}^(-1) - G_ss^(-1);
end

%% Plot poles
if 1

fn = fn + 1;
figure(fn)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.6]);

f1 = subplot(2,1,1);
for i = 1:length(pole_sys)
    scatter(real(pole_sys{i}),imag(pole_sys{i}),'x','LineWidth',1.5); hold on; grid on;
end
xlabel({'Real Part (Hz)'})
ylabel({'Imaginary Part (Hz)'})
X_H = 10;
X_L = -30;
XTick = [-30,-20,-10,0,10];
set(gca,'XLim',[X_L,X_H]);
set(gca,'XTick',XTick);
Y_H = 80;
Y_L = -Y_H;
YTick = [-80,-40,0,40,80];
set(gca,'YLim',[Y_L,Y_H]);
set(gca,'YTick',YTick);

f2 = subplot(2,1,2);
for i = 1:length(pole_sys)
    scatter(real(pole_sys{i}),imag(pole_sys{i}),'x','LineWidth',1.5); hold on; grid on;
end
xlabel('Real Part (Hz)')
ylabel({'Imaginary Part (Hz)'})
X_H = 0.006;
X_L = -0.006;
XTick = [-0.006,-0.003,0,0.003,0.006];
set(gca,'XLim',[X_L,X_H]);
set(gca,'XTick',XTick);
Y_H = 2;
Y_L = -Y_H;
YTick = [-2,-1,0,1,2];
set(gca,'YLim',[Y_L,Y_H]);
set(gca,'YTick',YTick);

if enable_save
    print(gcf,'Case_SG_PoleLocus.png','-dpng','-r600');
end

end

%% Plot K in complex plane
if 1

fn = fn+1;
figure(fn)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.25 0.35]);

plot_K(K{1},1j*2*pi); grid on; hold on; %
plot_K(K{2},1j*2*pi); grid on; hold on; % 
plot_K(K{3},1j*2*pi); grid on; hold on; % 
plot_K(K{4},1j*2*pi); grid on; hold on; % 
plot_K(K{5},1j*2*pi); grid on; hold on; %
X_L = -3e-6;
X_H = 2e-6;
set(gca,'XLim',[X_L,X_H]);
% set(gca,'XTick',XTick);
Y_L = -0.6e-3;
Y_H = 0.2e-3;
set(gca,'YLim',[Y_L,Y_H]);
set(gca,'YTick',[-0.6e-3,-0.4e-3,-0.2e-3,0,0.2e-3]);

if enable_save
    print(gcf,'Case_SG_Complex_K.png','-dpng','-r600');
end

end

%% Plot port
if 0
    
fn = fn+1;
figure(fn)
bodec(K{1},1j*omega_p,2*pi);
bodec(K{2},1j*omega_p,2*pi);
bodec(K{3},1j*omega_p,2*pi);
bodec(K{4},1j*omega_p,2*pi);
bodec(K{5},1j*omega_p,2*pi);
legend('1','2');

if enable_save
print(gcf,'Case_SG_K.png','-dpng','-r600');
end

end

if 0

fn = fn+1;
figure(fn)
bodec(Gop{1},1j*omega_p,2*pi);
bodec(Gop{6},1j*omega_p,2*pi);
legend('1','2');

if enable_save
print(gcf,'Case_SG_Gop.png','-dpng','-r600');
end

fn = fn+1;
figure(fn)
bodec(Gop{1},1j*omega_p,2*pi,'Color',color1);
bodec(Gop{6},1j*omega_p,2*pi,'Color',color2,'PhaseShift',-2*pi);
bodec(K{1},1j*omega_p,2*pi,'Color',color1);
bodec(K{6},1j*omega_p,2*pi,'Color',color2);
legend({'$R$=0pu','$R$=0.3pu'},'interpreter','latex');

if enable_save
print(gcf,'Case_SG_Gop_K.png','-dpng','-r600');
end

end

%% Plot transfer functions
if 0
    
% New
fn = fn+1;
figure(fn)
bodec(G_vi_c{i}(1,1),1j*omega_pn,2*pi,'PhaseOn',0);
bodec(G_vw_c{i}(1.1),1j*omega_pn,2*pi,'PhaseOn',0);
bodec(G_Ti_c{i}(1,1),1j*omega_pn,2*pi,'PhaseOn',0);
bodec(G_Tw{i},1j*omega_pn,2*pi,'PhaseOn',0);
% Set legend
legend({'$G_{vi}^{\prime\prime\prime}$','$G_{v\omega}^{\prime\prime\prime}$','$G_{Ti}^{\prime\prime\prime}$','$G_{T\omega}^{\prime\prime\prime}$'},'interpreter','latex','location','northeast')
Case1_Bode_Setting();
if enable_save
print(gcf,'Case_SG_Bode_TF.png','-dpng','-r600');
end

end

