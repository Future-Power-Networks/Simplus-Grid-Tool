% Author(s): Yitong Li

%%
clear all
clc
close all

%%
enable_save = 1;

xtick1 = 0.2;
xtick2 = 0.26;

%%
% Load data
Data_GFL1 = load('Data_GFL1').Data_GFL1;
Data_GFL2 = load('Data_GFL2').Data_GFL2;
Data_AngleDiff = load('Data_AngleDiff').Data_AngleDiff;

% time
time = Data_GFL1.time;
time_shift = 1.8;
time = time - time_shift;

% Organize data
vdq1   = Data_GFL1.signals(1).values;
idq1   = Data_GFL1.signals(2).values;
w1     = Data_GFL1.signals(7).values;

vdq2   = Data_GFL2.signals(1).values;
idq2   = Data_GFL2.signals(2).values;
w2     = Data_GFL2.signals(7).values;

AngleDiff     = Data_AngleDiff.signals(1).values;

%%
% Plot

fn = 0;
LineWidth = 1;
XLim = [0,0.6];
XTicks = [0,0.2,0.4,0.6];

% ### Grid following inverter 1
fn = fn+1;
figure(fn)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.18 0.6]); % position of left-bottow cornor + length/depth of figure

subplot(3,1,1)
plot(time,w1,'LineWidth',LineWidth); grid on; hold on;
ylabel('$\omega_1$ (pu)','interpreter','latex')
xlim(XLim);
xticks(XTicks);
ylim([0.9,1.2]);
yticks([0.9,1,1.1,1.2]);
plot([xtick1,xtick1],[-5,5],'--k'); hold on; grid on;
plot([xtick2,xtick2],[-5,5],'--k'); hold on; grid on;

subplot(3,1,2)
plot(time,idq1,'LineWidth',LineWidth); grid on; hold on;
ylabel('$i_{dq1}$ (pu)','interpreter','latex')
xlim(XLim);
xticks(XTicks);
%ylim([-4,1]);
%yticks([-4,-3,-2,-1,0,1]);
ylim([-1,4]);
yticks([-1,0,1,2,3,4]);
plot([xtick1,xtick1],[-5,5],'--k'); hold on; grid on;
plot([xtick2,xtick2],[-5,5],'--k'); hold on; grid on;

subplot(3,1,3)
plot(time,vdq1,'LineWidth',LineWidth); grid on; hold on;
ylabel('$v_{dq1}$ (pu)','interpreter','latex')
xlim(XLim);
xticks(XTicks);
ylim([-0.5,1.5]);
yticks([-0.5,0,0.5,1,1.5]);
xlabel('Time (s)','interpreter','latex')
plot([xtick1,xtick1],[-5,5],'--k'); hold on; grid on;
plot([xtick2,xtick2],[-5,5],'--k'); hold on; grid on;

if enable_save; print(gcf,'Sim_TwoGfl_Gfl1.png','-dpng','-r600'); end

% ### Grid following inverter 2
fn = fn+1;
figure(fn)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.18 0.6]); % position of left-bottow cornor + length/depth of figure

subplot(3,1,1)
plot(time,w2,'LineWidth',LineWidth); grid on; hold on;
ylabel('$\omega_2$ (pu)','interpreter','latex')
xlim(XLim);
xticks(XTicks);
ylim([0.9,1.2]);
yticks([0.9,1,1.1,1.2]);
plot([xtick1,xtick1],[-5,5],'--k'); hold on; grid on;
plot([xtick2,xtick2],[-5,5],'--k'); hold on; grid on;

subplot(3,1,2)
plot(time,idq2,'LineWidth',LineWidth); grid on; hold on;
ylabel('$i_{dq2}$ (pu)','interpreter','latex')
xlim(XLim);
xticks(XTicks);
% ylim([-4,1]);
% yticks([-4,-3,-2,-1,0,1]);
ylim([-1,4]);
yticks([-1,0,1,2,3,4]);
plot([xtick1,xtick1],[-5,5],'--k'); hold on; grid on;
plot([xtick2,xtick2],[-5,5],'--k'); hold on; grid on;

subplot(3,1,3)
plot(time,vdq2,'LineWidth',LineWidth); grid on; hold on;
ylabel('$v_{dq2}$ (pu)','interpreter','latex')
xlim(XLim);
xticks(XTicks);
ylim([-0.5,1.5]);
yticks([-0.5,0,0.5,1,1.5]);
xlabel('Time (s)','interpreter','latex')
plot([xtick1,xtick1],[-5,5],'--k'); hold on; grid on;
plot([xtick2,xtick2],[-5,5],'--k'); hold on; grid on;

if enable_save; print(gcf,'Sim_TwoGfl_Gfl2.png','-dpng','-r600'); end

% ### Angle difference
fn = fn+1;
figure(fn)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.26]); % position of left-bottow cornor + length/depth of figure

plot(time,AngleDiff,'LineWidth',LineWidth); grid on; hold on;
ylabel('$\theta_\Delta$ (degree)','interpreter','latex')
xlim(XLim);
xticks(XTicks);
ylim([-10,50]);
yticks([-10,10,30,50]);
plot([xtick1,xtick1],[-20,50],'--k'); hold on; grid on;
plot([xtick2,xtick2],[-20,50],'--k'); hold on; grid on;
xlabel('Time (s)','interpreter','latex')

if enable_save; print(gcf,'Sim_TwoGfl_AngleDiff.png','-dpng','-r600'); end
