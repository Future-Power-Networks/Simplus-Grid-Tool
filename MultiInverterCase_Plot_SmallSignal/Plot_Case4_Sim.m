% This function plots the simulation results of the single-SG-infinite-bus
% system.

% Author(s): Yitong Li

%%
close all
clear all
clc

%%
% Load simulation results
load('Case4_IBR1.mat');
load('Case4_IBR2.mat');
load('Case4_IBR3.mat');
load('Case4_IBR6.mat');
load('Case4_IBR8.mat');

% Organize data
fn = 0;
t_shift = -0.6;
t_vec = Case4_IBR1.time;
t_vec = t_vec + t_shift;

% Get w
w1_vec_ = Case4_IBR1.signals(5).values;
w2_vec_ = Case4_IBR2.signals(5).values;
w3_vec_ = Case4_IBR3.signals(5).values;
w6_vec_ = Case4_IBR6.signals(5).values;
w8_vec_ = Case4_IBR8.signals(5).values;

vdq1 = Case4_IBR1.signals(1).values;
vdq2 = Case4_IBR2.signals(1).values;
vdq3 = Case4_IBR3.signals(1).values;
vdq6 = Case4_IBR6.signals(1).values;
vdq8 = Case4_IBR8.signals(1).values;

idq1 = Case4_IBR1.signals(2).values;
idq2 = Case4_IBR2.signals(2).values;
idq3 = Case4_IBR3.signals(2).values;
idq6 = Case4_IBR6.signals(2).values;
idq8 = Case4_IBR8.signals(2).values;

N = length(w1_vec_);

for i = 1:N
    w1_vec(i) = w1_vec_(1,1,i);
    w2_vec(i) = w2_vec_(1,1,i);
    w3_vec(i) = w3_vec_(1,1,i);
    w6_vec(i) = w6_vec_(1,1,i);
    w8_vec(i) = w8_vec_(1,1,i);
end

clear('w1_vec_','w2_vec_','w3_vec_','w6_vec_','w8_vec_')

% Get pq
pq1_vec = Case4_IBR1.signals(7).values;
pq2_vec = Case4_IBR2.signals(7).values;
pq3_vec = Case4_IBR3.signals(7).values;
pq6_vec = Case4_IBR6.signals(7).values;
pq8_vec = Case4_IBR8.signals(7).values;

%% Plot
% Plot
fn = fn+1;
figure(fn)

rmax = 5;
cmax = 3;

% Set common
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.42 0.57]); % position of left-bottow cornor + length/depth of figure

X_L = 0;
X_H = 1.2;
XLim = [X_L,X_H];
XTick = [X_L,0.4,0.8,X_H];

Y_L_w = 0.85;
Y_H_w = 1.15;
YLim_w = [Y_L_w,Y_H_w];
YTick_w = [Y_L_w,1,Y_H_w];

YLim_i_GFM = [-8,8];
YTick_i_GFM = [-8,0,8];

YLim_i_GFL = [-1,1];
YTick_i_GFL = [-1,0,1];

Y_L_v = -0.2;
Y_H_v = 1.2;
YLim_v = [Y_L_v,Y_H_v];
YTick_v = [0,0.5,1];

xtick1 = 0.4;
xtick2 = 0.8;

Nw = -2;
Ni = -1;
Nv = 0;

LineWidth = 0.7;

% ### w
subplot(rmax,cmax,1*cmax+Nw);
plot(t_vec,w1_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',YLim_w);
set(gca,'YTick',YTick_w);
ylabel({'$\omega_1$ (pu)'},'interpreter','latex');
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

subplot(rmax,cmax,2*cmax+Nw);
plot(t_vec,w2_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',YLim_w);
set(gca,'YTick',YTick_w);
ylabel({'$\omega_2$ (pu)'},'interpreter','latex')
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

subplot(rmax,cmax,3*cmax+Nw);
plot(t_vec,w3_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',YLim_w);
set(gca,'YTick',YTick_w);
ylabel({'$\omega_3$ (pu)'},'interpreter','latex')
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

subplot(rmax,cmax,4*cmax+Nw);
plot(t_vec,w6_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',YLim_w);
set(gca,'YTick',YTick_w);
ylabel({'$\omega_6$ (pu)'},'interpreter','latex')
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

subplot(rmax,cmax,5*cmax+Nw);
plot(t_vec,w8_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',YLim_w);
set(gca,'YTick',YTick_w);
ylabel({'$\omega_8$ (pu)'},'interpreter','latex')
xlabel({'Time (s)'},'interpreter','latex');
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

% ### current
subplot(rmax,cmax,1*cmax+Ni);
plot(t_vec,idq1,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',YLim_i_GFM);
set(gca,'YTick',YTick_i_GFM);
ylabel({'$i_{dq1}$ (pu)'},'interpreter','latex');
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

subplot(rmax,cmax,2*cmax+Ni);
plot(t_vec,idq2,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',YLim_i_GFL);
set(gca,'YTick',YTick_i_GFL);
ylabel({'$i_{dq2}$ (pu)'},'interpreter','latex')
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

subplot(rmax,cmax,3*cmax+Ni);
plot(t_vec,idq3,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',YLim_i_GFM);
set(gca,'YTick',YTick_i_GFM);
ylabel({'$i_{dq3}$ (pu)'},'interpreter','latex')
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

subplot(rmax,cmax,4*cmax+Ni);
plot(t_vec,idq6,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',YLim_i_GFM);
set(gca,'YTick',YTick_i_GFM);
ylabel({'$i_{dq6}$ (pu)'},'interpreter','latex')
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

subplot(rmax,cmax,5*cmax+Ni);
plot(t_vec,idq8,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',YLim_i_GFL);
set(gca,'YTick',YTick_i_GFL);
ylabel({'$i_{dq8}$ (pu)'},'interpreter','latex')
xlabel({'Time (s)'},'interpreter','latex');
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

% ### voltage
subplot(rmax,cmax,1*cmax+Nv);
plot(t_vec,vdq1,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',YLim_v);
set(gca,'YTick',YTick_v);
ylabel({'$v_{dq1}$ (pu)'},'interpreter','latex');
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

subplot(rmax,cmax,2*cmax+Nv);
plot(t_vec,vdq2,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',YLim_v);
set(gca,'YTick',YTick_v);
ylabel({'$v_{dq2}$ (pu)'},'interpreter','latex')
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

subplot(rmax,cmax,3*cmax+Nv);
plot(t_vec,vdq3,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',YLim_v);
set(gca,'YTick',YTick_v);
ylabel({'$v_{dq3}$ (pu)'},'interpreter','latex')
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

subplot(rmax,cmax,4*cmax+Nv);
plot(t_vec,vdq6,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',YLim_v);
set(gca,'YTick',YTick_v);
ylabel({'$v_{dq6}$ (pu)'},'interpreter','latex')
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

subplot(rmax,cmax,5*cmax+Nv);
plot(t_vec,vdq8,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',YLim_v);
set(gca,'YTick',YTick_v);
ylabel({'$v_{dq8}$ (pu)'},'interpreter','latex')
xlabel({'Time (s)'},'interpreter','latex');
plot([xtick1,xtick1],YLim_w,'--k'); hold on; grid on;
plot([xtick2,xtick2],YLim_w,'--k'); hold on; grid on;

% f5 = subplot(3,2,5);
% f5.Position(1) = f5.Position(1) + 0.2;
% f5.Position(2) = f5.Position(2) + 0.02; 
% f5.Position(4) = f5.Position(4) - 0.015; 

if 1
    print(gcf,'Case4_GFM_Sim.png','-dpng','-r600');
end