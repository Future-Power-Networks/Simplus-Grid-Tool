% Author(s): Yitong Li

%%
clear all
clc
close all

%%
% Load data
Case3_PLL = load('Case3_PLL').Case3_PLL;

% time
time = Case3_PLL.time;
time_shift = 8.6;
time = time - time_shift;

% Organize data
vcdq_PLL = Case3_PLL.signals(1).values;
ildq_PLL = Case3_PLL.signals(2).values;
igdq_PLL = Case3_PLL.signals(3).values;
vcabc_PLL = Case3_PLL.signals(4).values;
ilabc_PLL = Case3_PLL.signals(5).values;
igabc_PLL = Case3_PLL.signals(6).values;
w_PLL = Case3_PLL.signals(7).values;
theta_PLL = Case3_PLL.signals(8).values;


%%
% Plot

fn = 0;
LineWidth = 1;
XLim = [0,1.2];
XTicks = [0,0.4,0.8,1.2];

YLim_vi = [-0.6,0.6];
YTicks_vi = [-0.5,0,0.5];

YLim_w = [0.8,1.4];
YTicks_w = [0.8,1,1.4];

XTick1 = [0.4,0.4];
XTick2 = [0.8,0.8];

fn = fn+1;
figure(fn)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.2 0.4]); % position of left-bottow cornor + length/depth of figure

subplot(3,1,1)
plot(time,vcdq_PLL,'LineWidth',LineWidth); hold on; grid on;
% plot(time,vcabc_PLL,'LineWidth',LineWidth); hold on; grid on;
ylabel('$v_{dq}$ (pu)','interpreter','latex')
ylim(YLim_vi);
yticks(YTicks_vi);
xlim(XLim);
xticks(XTicks);
% xlabel('Time (s)','interpreter','latex')
plot(XTick1,YLim_vi,'--k'); hold on; grid on;
plot(XTick2,YLim_vi,'--k'); hold on; grid on;

subplot(3,1,2)
plot(time,ildq_PLL,'LineWidth',LineWidth); hold on; grid on;
% plot(time,ilabc_PLL,'LineWidth',LineWidth); hold on; grid on;
ylabel('$i_{dq}$ (pu)','interpreter','latex')
ylim(YLim_vi);
yticks(YTicks_vi);
xlim(XLim);
xticks(XTicks);
% xlabel('Time (s)','interpreter','latex')
plot(XTick1,YLim_vi,'--k'); hold on; grid on;
plot(XTick2,YLim_vi,'--k'); hold on; grid on;

subplot(3,1,3)
plot(time,w_PLL,'LineWidth',LineWidth); hold on; grid on;
ylabel('$\omega$ (pu)','interpreter','latex')
xlim(XLim);
xticks(XTicks);
ylim(YLim_w);
yticks(YTicks_w);
xlabel('Time (s)','interpreter','latex')
plot(XTick1,YLim_w,'--k'); hold on; grid on;
plot(XTick2,YLim_w,'--k'); hold on; grid on;

if 0
    print(gcf,'Case3_Sim_GFL.png','-dpng','-r600'); 
end
