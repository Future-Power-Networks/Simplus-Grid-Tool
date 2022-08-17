% Author(s): Yitong Li

%%
clear all
clc
close all

%%
enable_save = 1;

xtick1 = 0.2;
xtick2 = 0.4;

%%
% Load data
Case2_GFM_Z = load('Case2_GFM_Z').Case2_GFM_Z;          % Grid-forming
Case2_GFL_Y = load('Case2_GFL_Y').Case2_GFL_Y;          % Grid-following

Case2_GFM_m     = load('Case2_GFM_m').Case2_GFM_m;          % Grid-forming
Case2_GFL_w_pll = load('Case2_GFL_w_pll').Case2_GFL_w_pll;	% Grid-following

Case2_GFM_w_v = load('Case2_GFM_w_v').Case2_GFM_w_v;    % Grid-forming
Case2_GFL_w_i = load('Case2_GFL_w_i').Case2_GFL_w_i;    % Grid-following


% time
time = Case2_GFM_Z.time;
time_shift = 1.8;
time = time - time_shift;

% Organize data
w_GFM_Z     = Case2_GFM_Z.signals(7).values;
w_GFL_Y     = Case2_GFL_Y.signals(7).values;

w_GFM_m     = Case2_GFM_m.signals(7).values;
w_GFL_w_pll = Case2_GFL_w_pll.signals(7).values;

w_GFM_w_v   = Case2_GFM_w_v.signals(7).values;
w_GFL_w_i   = Case2_GFL_w_i.signals(7).values;

%%
% Plot

fn = 0;
LineWidth = 1;
XLim = [0,0.6];
XTicks = [0,0.2,0.4,0.6];

% ### Grid forming
% Change Z
fn = fn+1;
figure(fn)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.2 0.25]); % position of left-bottow cornor + length/depth of figure
plot(time,w_GFM_Z,'LineWidth',LineWidth); grid on; hold on;
ylabel('$\omega$ (pu)','interpreter','latex')
xlim(XLim);
xticks(XTicks);
ylim([0.999,1.001]);
yticks([0.999,1,1.001]);
xlabel('Time (s)','interpreter','latex')
plot([xtick1,xtick1],[0.999,1.001],'--k'); hold on; grid on;
plot([xtick2,xtick2],[0.999,1.001],'--k'); hold on; grid on;
if enable_save; print(gcf,'Case2_Sim_GFM_Z.png','-dpng','-r600'); end

% Change m
fn = fn+1;
figure(fn)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.2 0.25]); % position of left-bottow cornor + length/depth of figure
plot(time,w_GFM_m,'LineWidth',LineWidth); grid on; hold on;
ylabel('$\omega$ (pu)','interpreter','latex')
xlim(XLim);
xticks(XTicks);
ylim([0.96,1.04]);
yticks([0.96,1,1.04]);
xlabel('Time (s)','interpreter','latex')
plot([xtick1,xtick1],[0.96,1.04],'--k'); hold on; grid on;
plot([xtick2,xtick2],[0.96,1.04],'--k'); hold on; grid on;
if enable_save; print(gcf,'Case2_Sim_GFM_m.png','-dpng','-r600'); end

% Change w_v
fn = fn+1;
figure(fn)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.2 0.25]); % position of left-bottow cornor + length/depth of figure
plot(time,w_GFM_w_v,'LineWidth',LineWidth); grid on; hold on;
ylabel('$\omega$ (pu)','interpreter','latex')
xlim(XLim);
xticks(XTicks);
ylim([0.9999,1.0001]);
yticks([0.9999,1,1.0001]);
xlabel('Time (s)','interpreter','latex')
plot([xtick1,xtick1],[0.9999,1.0001],'--k'); hold on; grid on;
plot([xtick2,xtick2],[0.9999,1.0001],'--k'); hold on; grid on;
if enable_save; print(gcf,'Case2_Sim_GFM_w_v.png','-dpng','-r600'); end

% ### Grid following 
% Change Y
fn = fn+1;
figure(fn)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.2 0.25]); % position of left-bottow cornor + length/depth of figure
plot(time,w_GFL_Y,'LineWidth',LineWidth); grid on; hold on;
ylabel('$\omega$ (pu)','interpreter','latex')
xlim(XLim);
xticks(XTicks);
ylim([0.85,1.15]);
yticks([0.85,1,1.15]);
xlabel('Time (s)','interpreter','latex')
plot([xtick1,xtick1],[0.85,1.15],'--k'); hold on; grid on;
plot([xtick2,xtick2],[0.85,1.15],'--k'); hold on; grid on;
if enable_save; print(gcf,'Case2_Sim_GFL_Y.png','-dpng','-r600');  end

% Change w_pll
fn = fn+1;
figure(fn)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.2 0.25]); % position of left-bottow cornor + length/depth of figure
plot(time,w_GFL_w_pll,'LineWidth',LineWidth); grid on; hold on;
ylabel('$\omega$ (pu)','interpreter','latex')
xlim(XLim);
xticks(XTicks);
ylim([0.98,1.02]);
yticks([0.98,1,1.02]);
xlabel('Time (s)','interpreter','latex')
plot([xtick1,xtick1],[0.98,1.02],'--k'); hold on; grid on;
plot([xtick2,xtick2],[0.98,1.02],'--k'); hold on; grid on;
if enable_save; print(gcf,'Case2_Sim_GFL_w_pll.png','-dpng','-r600');  end

% Change w_i
fn = fn+1;
figure(fn)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.2 0.25]); % position of left-bottow cornor + length/depth of figure
plot(time,w_GFL_w_i,'LineWidth',LineWidth); grid on; hold on;
ylabel('$\omega$ (pu)','interpreter','latex')
xlim(XLim);
xticks(XTicks);
ylim([0.9998,1.0002]);
yticks([0.9998,1,1.0002]);
xlabel('Time (s)','interpreter','latex')
plot([xtick1,xtick1],[0.9998,1.0002],'--k'); hold on; grid on;
plot([xtick2,xtick2],[0.9998,1.0002],'--k'); hold on; grid on;
if enable_save; print(gcf,'Case2_Sim_GFL_w_i.png','-dpng','-r600');  end
