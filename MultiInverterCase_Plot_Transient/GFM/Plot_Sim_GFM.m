% This function plots the simulation results of the single-SG-infinite-bus
% system.

% Author(s): Yitong Li

%%
close all
clear all
clc

RgbBlue = [0, 0.4470, 0.7410];

%%
% Load simulation results
Stable_Data_App1 = load('Stable_Data_App1.mat').Data_App1;
Stable_Data_App2 = load('Stable_Data_App2.mat').Data_App2;
Stable_Data_App3 = load('Stable_Data_App3.mat').Data_App3;
Stable_Data_App6 = load('Stable_Data_App6.mat').Data_App6;
Stable_Data_App8 = load('Stable_Data_App8.mat').Data_App8;

Unstable_Data_App1 = load('Unstable_Data_App1.mat').Data_App1;
Unstable_Data_App2 = load('Unstable_Data_App2.mat').Data_App2;
Unstable_Data_App3 = load('Unstable_Data_App3.mat').Data_App3;
Unstable_Data_App6 = load('Unstable_Data_App6.mat').Data_App6;
Unstable_Data_App8 = load('Unstable_Data_App8.mat').Data_App8;
% Unstable_Data_VoltageBus6 = load('Unstable_Data_VoltageBus6.mat').Data_VoltageBus6;

Load_Data_App1 = load('Load_Data_App1.mat').Data_App1;
Load_Data_App2 = load('Load_Data_App2.mat').Data_App2;
Load_Data_App3 = load('Load_Data_App3.mat').Data_App3;
Load_Data_App6 = load('Load_Data_App6.mat').Data_App6;
Load_Data_App8 = load('Load_Data_App8.mat').Data_App8;

% Organize data
fn = 0;
t_shift = -4.8;
t_vec = Stable_Data_App1.time;
t_vec = t_vec + t_shift;

% % Get voltage
% st_vabc = Stable_Data_VoltageBus6.signals(1).values;
% st_vdq = Stable_Data_VoltageBus6.signals(2).values;
% 
% us_vabc = Unstable_Data_VoltageBus6.signals(1).values;
% us_vdq = Unstable_Data_VoltageBus6.signals(2).values;

% Get w
w1_vec_ = Stable_Data_App1.signals(5).values;
w2_vec_ = Stable_Data_App2.signals(5).values;
w3_vec_ = Stable_Data_App3.signals(5).values;
w6_vec_ = Stable_Data_App6.signals(5).values;
w8_vec_ = Stable_Data_App8.signals(5).values;
N = length(w1_vec_);
for i = 1:N
    st_w1_vec(i) = w1_vec_(1,1,i);
    st_w2_vec(i) = w2_vec_(1,1,i);
    st_w3_vec(i) = w3_vec_(1,1,i);
    st_w6_vec(i) = w6_vec_(1,1,i);
    st_w8_vec(i) = w8_vec_(1,1,i);
end
clear('w1_vec_','w2_vec_','w3_vec_','w6_vec_','w8_vec_');

w1_vec_ = Unstable_Data_App1.signals(5).values;
w2_vec_ = Unstable_Data_App2.signals(5).values;
w3_vec_ = Unstable_Data_App3.signals(5).values;
w6_vec_ = Unstable_Data_App6.signals(5).values;
w8_vec_ = Unstable_Data_App8.signals(5).values;
N = length(w1_vec_);
for i = 1:N
    us_w1_vec(i) = w1_vec_(1,1,i);
    us_w2_vec(i) = w2_vec_(1,1,i);
    us_w3_vec(i) = w3_vec_(1,1,i);
    us_w6_vec(i) = w6_vec_(1,1,i);
    us_w8_vec(i) = w8_vec_(1,1,i);
end
clear('w1_vec_','w2_vec_','w3_vec_','w6_vec_','w8_vec_');

w1_vec_ = Load_Data_App1.signals(5).values;
w2_vec_ = Load_Data_App2.signals(5).values;
w3_vec_ = Load_Data_App3.signals(5).values;
w6_vec_ = Load_Data_App6.signals(5).values;
w8_vec_ = Load_Data_App8.signals(5).values;
N = length(w1_vec_);
for i = 1:N
    ld_w1_vec(i) = w1_vec_(1,1,i);
    ld_w2_vec(i) = w2_vec_(1,1,i);
    ld_w3_vec(i) = w3_vec_(1,1,i);
    ld_w6_vec(i) = w6_vec_(1,1,i);
    ld_w8_vec(i) = w8_vec_(1,1,i);
end
clear('w1_vec_','w2_vec_','w3_vec_','w6_vec_','w8_vec_');

%% Plot
% Plot
fn = fn+1;
figure(fn)

rmax = 5;
cmax = 1;

% Set common
% [0.1 0.1 0.42 0.57]


X_L = 0.1;
X_H = 0.6;
XLim = [X_L,X_H];
xtick1 = 0.2;
xtick2 = 0.26;
XTick = [X_L,xtick1,X_H];

Y_L_w = -2;
Y_H_w = 2;
YLim_w = [Y_L_w,Y_H_w];

LineWidth = 1;

% % ### voltage
% r=r+1;
% subplot(rmax,cmax,r);
% plot(t_vec,st_vdq(:,1),'LineWidth',LineWidth); grid on; hold on;
% plot(t_vec,us_vdq(:,1),'-.','LineWidth',LineWidth,'Color',RgbBlue); grid on; hold on;
% set(gca,'XLim',XLim);
% set(gca,'XTick',XTick);
% set(gca,'YLim',[0,1.02]);
% set(gca,'YTick',[0.8,1,1.02]);
% ylabel({'$v_{d6}$ (pu)'},'interpreter','latex');
% plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;

%%
figure(1)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.2 0.7]); % position of left-bottow cornor + length/depth of figure
r = 0;

% ### w
r=r+1;
subplot(rmax,cmax,r);
plot(t_vec,st_w1_vec,'-.','LineWidth',LineWidth); grid on; hold on;
plot(t_vec,us_w1_vec,'LineWidth',LineWidth,'Color',RgbBlue); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',[0.9,1.1]);
set(gca,'YTick',[0.9,1,1.1]);
ylabel({'$\omega_1$ (pu)'},'interpreter','latex');
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

r=r+1;
subplot(rmax,cmax,r);
plot(t_vec,st_w2_vec,'-.','LineWidth',LineWidth); grid on; hold on;
plot(t_vec,us_w2_vec,'LineWidth',LineWidth,'Color',RgbBlue); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',[0.9,1.1]);
set(gca,'YTick',[0.9,1,1.1]);
ylabel({'$\omega_2$ (pu)'},'interpreter','latex')
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

r=r+1;
subplot(rmax,cmax,r);
plot(t_vec,st_w3_vec,'-.','LineWidth',LineWidth); grid on; hold on;
plot(t_vec,us_w3_vec,'LineWidth',LineWidth,'Color',RgbBlue); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',[0.9,1.1]);
set(gca,'YTick',[0.9,1,1.1]);
ylabel({'$\omega_3$ (pu)'},'interpreter','latex')
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

r=r+1;
subplot(rmax,cmax,r);
plot(t_vec,st_w6_vec,'-.','LineWidth',LineWidth); grid on; hold on;
plot(t_vec,us_w6_vec,'LineWidth',LineWidth,'Color',RgbBlue); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',[-1.3,1.1]);
set(gca,'YTick',[-1.3,1.1]);
ylabel({'$\omega_6$ (pu)'},'interpreter','latex')
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

r=r+1;
subplot(rmax,cmax,r);
plot(t_vec,st_w8_vec,'-.','LineWidth',LineWidth); grid on; hold on;
plot(t_vec,us_w8_vec,'LineWidth',LineWidth,'Color',RgbBlue); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',[0.9,1.1]);
set(gca,'YTick',[0.9,1,1.1]);
ylabel({'$\omega_8$ (pu)'},'interpreter','latex')
xlabel({'Time (s)'},'interpreter','latex');
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

if 1
    print(gcf,'Sim_Transient_GFM.png','-dpng','-r600');
end

%%
figure(2)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.2 0.7]); % position of left-bottow cornor + length/depth of figure
r = 0;

% ### w
r=r+1;
subplot(rmax,cmax,r);
plot(t_vec,ld_w1_vec,'-.','LineWidth',LineWidth); grid on; hold on;
plot(t_vec,us_w1_vec,'LineWidth',LineWidth,'Color',RgbBlue); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',[0.9,1.1]);
set(gca,'YTick',[0.9,1,1.1]);
ylabel({'$\omega_1$ (pu)'},'interpreter','latex');
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

r=r+1;
subplot(rmax,cmax,r);
plot(t_vec,ld_w2_vec,'-.','LineWidth',LineWidth); grid on; hold on;
plot(t_vec,us_w2_vec,'LineWidth',LineWidth,'Color',RgbBlue); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',[0.9,1.1]);
set(gca,'YTick',[0.9,1,1.1]);
ylabel({'$\omega_2$ (pu)'},'interpreter','latex')
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

r=r+1;
subplot(rmax,cmax,r);
plot(t_vec,ld_w3_vec,'-.','LineWidth',LineWidth); grid on; hold on;
plot(t_vec,us_w3_vec,'LineWidth',LineWidth,'Color',RgbBlue); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',[0.9,1.1]);
set(gca,'YTick',[0.9,1,1.1]);
ylabel({'$\omega_3$ (pu)'},'interpreter','latex')
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

r=r+1;
subplot(rmax,cmax,r);
plot(t_vec,ld_w6_vec,'-.','LineWidth',LineWidth); grid on; hold on;
plot(t_vec,us_w6_vec,'LineWidth',LineWidth,'Color',RgbBlue); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',[-1.3,1.1]);
set(gca,'YTick',[-1.3,1.1]);
ylabel({'$\omega_6$ (pu)'},'interpreter','latex')
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

r=r+1;
subplot(rmax,cmax,r);
plot(t_vec,ld_w8_vec,'-.','LineWidth',LineWidth); grid on; hold on;
plot(t_vec,us_w8_vec,'LineWidth',LineWidth,'Color',RgbBlue); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',[0.9,1.1]);
set(gca,'YTick',[0.9,1,1.1]);
ylabel({'$\omega_8$ (pu)'},'interpreter','latex')
xlabel({'Time (s)'},'interpreter','latex');
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

if 1
    print(gcf,'Sim_Transient_GFM_Load.png','-dpng','-r600');
end