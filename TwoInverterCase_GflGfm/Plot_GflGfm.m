% Author(s): Yitong Li

%%
clear all
clc
close all

%%
enable_save = 1;

%%
% Load data
Data_GFL1 = load('Data_GFL1').Data_GFL1;
Data_GFM2 = load('Data_GFM2').Data_GFM2;
Data_AngleDiff = load('Data_AngleDiff').Data_AngleDiff;

% time
time = Data_GFL1.time;
time_shift = 4;
time = time - time_shift;

% Organize data
vdq1   = Data_GFL1.signals(1).values;
idq1   = Data_GFL1.signals(2).values;
w1     = Data_GFL1.signals(7).values;

vdq2   = Data_GFM2.signals(1).values;
idq2   = Data_GFM2.signals(2).values;
w2     = Data_GFM2.signals(7).values;

AngleDiff = Data_AngleDiff.signals(1).values;

%%
% Plot

fn = 0;
LineWidth = 1;
XLim = [0,5];
XTicks = [0,1,2,3,4,5];

fn = fn+1;
figure(fn)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.28 0.55]); % position of left-bottow cornor + length/depth of figure

subplot(3,1,1)
plot(time,w1,'LineWidth',LineWidth); grid on; hold on;
plot([1,1],[1,1.1],'--k'); hold on; grid on;
ylabel('$\omega_1$ (pu)','interpreter','latex')
xlim(XLim);
xticks(XTicks);
ylim([1.04,1.06]);
yticks([1.04,1.05,1.06]);

subplot(3,1,2)
plot(time,w2,'LineWidth',LineWidth); grid on; hold on;
plot([1,1],[1,1.1],'--k'); hold on; grid on;
ylabel('$\omega_2$ (pu)','interpreter','latex')
xlim(XLim);
xticks(XTicks);
ylim([1.04,1.06]);
yticks([1.04,1.05,1.06]);

subplot(3,1,3)
plot(time,AngleDiff,'LineWidth',LineWidth); grid on; hold on;
plot([1,1],[-30,30],'--k'); hold on; grid on;
ylabel('$\theta_\Delta$ (degree)','interpreter','latex')
xlim(XLim);
xticks(XTicks);
ylim([-30,30]);
yticks([-30,-15,0,15,30]);
xlabel('Time (s)','interpreter','latex')

if enable_save; print(gcf,'Sim_GflGfm_Raw.png','-dpng','-r600'); end


