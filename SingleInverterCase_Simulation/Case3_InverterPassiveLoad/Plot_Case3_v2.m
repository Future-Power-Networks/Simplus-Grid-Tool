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
LineWidth = 0.7;
XLim = [0,1.2];
XTicks = [0,0.4,0.8,1.2];

YLim_v = [-0.2,0.7];
YTicks_v = [0,0.5];

YLim_i = [-0.7,0.2];
YTicks_i = [-0.5,0];

YLim_vq = [-0.1,0.2];
YTicks_vq = [-0.1,0,0.1,0.2];

YLim_iq = [-0.2,0.2];
YTicks_iq = [-0.2,0,0.2];

YLim_w = [0.8,1.4];
YTicks_w = [0.8,1,1.4];

XTick1 = [0.4,0.4];
XTick2 = [0.8,0.8];

DefaultRed = [0.8500, 0.3250, 0.0980];

fn = fn+1;
figure(fn)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]); % position of left-bottow cornor + length/depth of figure

subplot(3,2,1)
plot(time,w_PLL,'LineWidth',LineWidth); hold on; grid on;
ylabel('$\omega$ (pu)','interpreter','latex')
xlim(XLim);
xticks(XTicks);
ylim(YLim_w);
yticks(YTicks_w);
plot(XTick1,YLim_w,'--k'); hold on; grid on;
plot(XTick2,YLim_w,'--k'); hold on; grid on;

subplot(3,2,3)
plot(time,ildq_PLL,'LineWidth',LineWidth); hold on; grid on;
% plot(time,ilabc_PLL,'LineWidth',LineWidth); hold on; grid on;
ylabel('$i_{dq}$ (pu)','interpreter','latex')
ylim(YLim_i);
yticks(YTicks_i);
xlim(XLim);
xticks(XTicks);
% xlabel('Time (s)','interpreter','latex')
plot(XTick1,YLim_i,'--k'); hold on; grid on;
plot(XTick2,YLim_i,'--k'); hold on; grid on;

subplot(3,2,4)
plot(time,ildq_PLL(:,2),'LineWidth',LineWidth,'Color',DefaultRed); hold on; grid on;
% plot(time,ilabc_PLL,'LineWidth',LineWidth); hold on; grid on;
% ylabel('$i_{q}$ (pu)','interpreter','latex')
ylim(YLim_iq);
yticks(YTicks_iq);
xlim(XLim);
xticks(XTicks);
% xlabel('Time (s)','interpreter','latex')
plot(XTick1,YLim_i,'--k'); hold on; grid on;
plot(XTick2,YLim_i,'--k'); hold on; grid on;

subplot(3,2,5)
plot(time,vcdq_PLL,'LineWidth',LineWidth); hold on; grid on;
% plot(time,vcabc_PLL,'LineWidth',LineWidth); hold on; grid on;
ylabel('$v_{dq}$ (pu)','interpreter','latex')
ylim(YLim_v);
yticks(YTicks_v);
xlim(XLim);
xticks(XTicks);
% xlabel('Time (s)','interpreter','latex')
plot(XTick1,YLim_v,'--k'); hold on; grid on;
plot(XTick2,YLim_v,'--k'); hold on; grid on;
xlabel('Time (s)','interpreter','latex')

subplot(3,2,6)
plot(time,vcdq_PLL(:,2),'LineWidth',LineWidth,'Color',DefaultRed); hold on; grid on;
% plot(time,vcabc_PLL,'LineWidth',LineWidth); hold on; grid on;
% ylabel('$v_{q}$ (pu)','interpreter','latex')
ylim(YLim_vq);
yticks(YTicks_vq);
xlim(XLim);
xticks(XTicks);
% xlabel('Time (s)','interpreter','latex')
plot(XTick1,YLim_v,'--k'); hold on; grid on;
plot(XTick2,YLim_v,'--k'); hold on; grid on;
xlabel('Time (s)','interpreter','latex')

if 1
    print(gcf,'Case3_Sim_GFL.png','-dpng','-r600'); 
end
