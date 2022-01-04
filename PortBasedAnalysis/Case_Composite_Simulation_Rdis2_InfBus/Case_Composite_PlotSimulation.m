% This function plots the simulation results of the single-SG-infinite-bus
% system.

% Author(s): Yitong Li

%%
close all
clear all
clc

enable_save = 1;
FigSize = [0.1 0.1 0.28 0.7];
LineWidth = 1;

%%
% Load simulation results
load('Case_Composite_InfBus_SM1.mat');
load('Case_Composite_InfBus_SM2.mat');
load('Case_Composite_InfBus_SM6.mat');
load('Case_Composite_InfBus_VSI3.mat');
load('Case_Composite_InfBus_VSI8.mat');

% Organize data
fn = 0;
t_shift = -15;
t_vec = Case_Composite_SM1.time;
t_vec = t_vec + t_shift;

% Get w
w1_vec_ = Case_Composite_SM1.signals(5).values;
w2_vec_ = Case_Composite_SM2.signals(5).values;
w6_vec_ = Case_Composite_SM6.signals(5).values;

w3_vec_ = Case_Composite_VSI3.signals(5).values;
w8_vec_ = Case_Composite_VSI8.signals(5).values;

N = length(w2_vec_);

w1_vec = 0;
w2_vec = 0;
w6_vec = 0;

Wbase = 2*pi*60;    
Wbase_ = Wbase*0.9994; 	% The actual base used in inf bus case, in 
                      	% order to match the frequency of this case to 
                       	% the SG case
for i = 1:N
    w1_vec(i) = w1_vec_(1,1,i)*Wbase_/Wbase;
    w2_vec(i) = w2_vec_(1,1,i)*Wbase_/Wbase;
    w6_vec(i) = w6_vec_(1,1,i)*Wbase_/Wbase;
    
    w3_vec(i) = w3_vec_(1,1,i)*Wbase_/Wbase;
    w8_vec(i) = w8_vec_(1,1,i)*Wbase_/Wbase;
end

% Get pq
pq1_vec = Case_Composite_SM1.signals(7).values;
pq2_vec = Case_Composite_SM2.signals(7).values;
pq6_vec = Case_Composite_SM6.signals(7).values;

pq3_vec = Case_Composite_VSI3.signals(7).values;
pq8_vec = Case_Composite_VSI8.signals(7).values;

% Get pqfil
pqfil1_vec = Case_Composite_SM1.signals(8).values;
pqfil2_vec = Case_Composite_SM2.signals(8).values;
pqfil6_vec = Case_Composite_SM6.signals(8).values;

pqfil3_vec = Case_Composite_VSI3.signals(8).values;
pqfil8_vec = Case_Composite_VSI8.signals(8).values;

%% Plot
if 0

% Plot
fn = fn+1;
figure(fn)

% Set common
set(gcf,'units','normalized','outerposition',FigSize); % position of left-bottow cornor + length/depth of figure

X_L = -1;
X_H = 10;
XLim = [X_L,X_H];
XTick = [-1,0,5,10];

rmax = 5;
cmax = 1;

% ### Subplot
subplot(rmax,cmax,1);
plot(t_vec,w1_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
dY = 0.003;
Y_M = 0.9994;
Y_L = Y_M - dY;
Y_H = Y_M + dY;
set(gca,'YLim',[Y_L,Y_H]);
set(gca,'YTick',[Y_L,Y_M,Y_H]);
ylabel({'$\omega$ (pu)'},'interpreter','latex');
title('SG1')

% ### Subplot
subplot(rmax,cmax,1+cmax);
plot(t_vec,w2_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
dY = 0.003;
Y_M = 0.9994;
Y_L = Y_M - dY;
Y_H = Y_M + dY;
set(gca,'YLim',[Y_L,Y_H]);
set(gca,'YTick',[Y_L,Y_M,Y_H]);
ylabel({'$\omega$ (pu)'},'interpreter','latex')
title('SG2')

% ### Subplot
subplot(rmax,cmax,1+cmax*2);
plot(t_vec,w6_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
dY = 0.003;
Y_M = 0.9994;
Y_L = Y_M - dY;
Y_H = Y_M + dY;
set(gca,'YLim',[Y_L,Y_H]);
set(gca,'YTick',[Y_L,Y_M,Y_H]);
ylabel({'$\omega$ (pu)'},'interpreter','latex')
title('SG8')

if rmax == 5
% ### Subplot
subplot(rmax,cmax,1+cmax*3);
plot(t_vec,w3_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
dY = 0.003;
Y_M = 0.9994;
Y_L = Y_M - dY;
Y_H = Y_M + dY;
set(gca,'YLim',[Y_L,Y_H]);
set(gca,'YTick',[Y_L,Y_M,Y_H]);
ylabel({'$\omega$ (pu)'},'interpreter','latex')
title('VSI3')

% ### Subplot
subplot(rmax,cmax,1+cmax*4);
plot(t_vec,w8_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
dY = 0.003;
Y_M = 0.9994;
Y_L = Y_M - dY;
Y_H = Y_M + dY;
set(gca,'YLim',[Y_L,Y_H]);
set(gca,'YTick',[Y_L,Y_M,Y_H]);
ylabel({'$\omega$ (pu)'},'interpreter','latex')
title('VSI8')
end

xlabel({'Time (s)'},'interpreter','latex');

% if cmax >= 2
% % ### Subplot
% subplot(rmax,cmax,2);
% plot(t_vec,pqfil1_vec(:,1),'LineWidth',LineWidth); grid on; hold on;
% set(gca,'XLim',XLim);
% set(gca,'XTick',XTick);
% ylabel({'$P$ (pu)'},'interpreter','latex');
% title('SG1')
% 
% % ### Subplot
% subplot(rmax,cmax,2+cmax);
% plot(t_vec,pqfil2_vec(:,1),'LineWidth',LineWidth); grid on; hold on;
% set(gca,'XLim',XLim);
% set(gca,'XTick',XTick);
% ylabel({'$P$ (pu)'},'interpreter','latex');
% title('SG1')
% 
% % ### Subplot
% subplot(rmax,cmax,2+cmax*2);
% plot(t_vec,pqfil3_vec(:,1),'LineWidth',LineWidth); grid on; hold on;
% set(gca,'XLim',XLim);
% set(gca,'XTick',XTick);
% ylabel({'$P$ (pu)'},'interpreter','latex');
% title('SG1')
% end

if enable_save; print(gcf,'Case_Composite_Simulation.png','-dpng','-r600'); end

end

%% Plot zoomed

if 0
% Plot
fn = fn+1;
figure(fn)

% Set common
set(gcf,'units','normalized','outerposition',FigSize); % position of left-bottow cornor + length/depth of figure

X_L = 4;
X_H = 6;
XLim = [X_L,X_H];
XTick = [4,5,6];

rmax = 5;
cmax = 1;

% ### Subplot
subplot(rmax,cmax,1);
plot(t_vec,w1_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
dY = 0.0003;
Y_M = 0.9994;
Y_L = Y_M - dY;
Y_H = Y_M + dY;
set(gca,'YLim',[Y_L,Y_H]);
set(gca,'YTick',[Y_L,Y_M,Y_H]);
%xlabel({'Time (s)'},'interpreter','latex');
ylabel({'$\omega$ (pu)'},'interpreter','latex');
title('SG1')

% ### Subplot
subplot(rmax,cmax,1+cmax);
plot(t_vec,w2_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
dY = 0.003;
Y_M = 0.9994;
Y_L = Y_M - dY;
Y_H = Y_M + dY;
set(gca,'YLim',[Y_L,Y_H]);
set(gca,'YTick',[Y_L,Y_M,Y_H]);
%xlabel({'Time (s)'},'interpreter','latex');
ylabel({'$\omega$ (pu)'},'interpreter','latex')
title('SG2')

% ### Subplot
subplot(rmax,cmax,1+cmax*2);
plot(t_vec,w6_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
dY = 0.0003;
Y_M = 0.9994;
Y_L = Y_M - dY;
Y_H = Y_M + dY;
set(gca,'YLim',[Y_L,Y_H]);
set(gca,'YTick',[Y_L,Y_M,Y_H]);
ylabel({'$\omega$ (pu)'},'interpreter','latex')
title('SG8')

if rmax == 5
% ### Subplot
subplot(rmax,cmax,1+cmax*3);
plot(t_vec,w3_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
dY = 0.003;
Y_M = 0.9994;
Y_L = Y_M - dY;
Y_H = Y_M + dY;
set(gca,'YLim',[Y_L,Y_H]);
set(gca,'YTick',[Y_L,Y_M,Y_H]);
ylabel({'$\omega$ (pu)'},'interpreter','latex')
title('VSI3')

% ### Subplot
subplot(rmax,cmax,1+cmax*4);
plot(t_vec,w8_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
dY = 0.003;
Y_M = 0.9994;
Y_L = Y_M - dY;
Y_H = Y_M + dY;
set(gca,'YLim',[Y_L,Y_H]);
set(gca,'YTick',[Y_L,Y_M,Y_H]);
ylabel({'$\omega$ (pu)'},'interpreter','latex')
title('VSI8')
end

xlabel({'Time (s)'},'interpreter','latex');

if enable_save; print(gcf,'Case_Composite_Simulation_Zoom.png','-dpng','-r600'); end

end

%% Plot both
if 1
% Plot
fn = fn+1;
figure(fn)

FigSize = [0.1 0.1 0.34 0.8];

% Set common
set(gcf,'units','normalized','outerposition',FigSize); % position of left-bottow cornor + length/depth of figure

X_L = -1;
X_H = 10;
XLim = [X_L,X_H];
XTick = [-1,0,5,10];

rmax = 5;
cmax = 2;

% ### Subplot
subplot(rmax,cmax,1);
plot(t_vec,w1_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
dY = 0.0015;
Y_M = 0.9994;
Y_L = Y_M - dY;
Y_H = Y_M + dY;
set(gca,'YLim',[Y_L,Y_H]);
set(gca,'YTick',[Y_L,Y_M,Y_H]);
ylabel({'$\omega_1$ (pu)'},'interpreter','latex');
title('SG1')

% ### Subplot
subplot(rmax,cmax,1+cmax);
plot(t_vec,w2_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
dY = 0.003;
Y_M = 0.9994;
Y_L = Y_M - dY;
Y_H = Y_M + dY;
set(gca,'YLim',[Y_L,Y_H]);
set(gca,'YTick',[Y_L,Y_M,Y_H]);
ylabel({'$\omega_2$ (pu)'},'interpreter','latex')
title('SG2')

% ### Subplot
subplot(rmax,cmax,1+cmax*2);
plot(t_vec,w6_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
dY = 0.0015;
Y_M = 0.9994;
Y_L = Y_M - dY;
Y_H = Y_M + dY;
set(gca,'YLim',[Y_L,Y_H]);
set(gca,'YTick',[Y_L,Y_M,Y_H]);
ylabel({'$\omega_6$ (pu)'},'interpreter','latex')
title('SG6')

% ### Subplot
subplot(rmax,cmax,1+cmax*3);
plot(t_vec,w3_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
dY = 0.003;
Y_M = 0.9994;
Y_L = Y_M - dY;
Y_H = Y_M + dY;
set(gca,'YLim',[Y_L,Y_H]);
set(gca,'YTick',[Y_L,Y_M,Y_H]);
ylabel({'$\omega_3$ (pu)'},'interpreter','latex')
title('VSI3')

% ### Subplot
subplot(rmax,cmax,1+cmax*4);
plot(t_vec,w8_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
dY = 0.003;
Y_M = 0.9994;
Y_L = Y_M - dY;
Y_H = Y_M + dY;
set(gca,'YLim',[Y_L,Y_H]);
set(gca,'YTick',[Y_L,Y_M,Y_H]);
ylabel({'$\omega_8$ (pu)'},'interpreter','latex')
title('VSI8')

xlabel({'Time (s)'},'interpreter','latex');

X_L = -0.1;
X_H = 0.1;
XLim = [X_L,X_H];
XTick = [-0.1,0,0.1];

% ### Subplot
subplot(rmax,cmax,2);
plot(t_vec,w1_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
dY = 0.0003;
Y_M = 0.9994;
Y_L = Y_M - dY;
Y_H = Y_M + dY;
set(gca,'YLim',[Y_L,Y_H]);
set(gca,'YTick',[Y_L,Y_M,Y_H]);

% ### Subplot
subplot(rmax,cmax,2+cmax);
plot(t_vec,w2_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
dY = 0.003;
Y_M = 0.9994;
Y_L = Y_M - dY;
Y_H = Y_M + dY;
set(gca,'YLim',[Y_L,Y_H]);
set(gca,'YTick',[Y_L,Y_M,Y_H]);

% ### Subplot
subplot(rmax,cmax,2+cmax*2);
plot(t_vec,w6_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
dY = 0.0003;
Y_M = 0.9994;
Y_L = Y_M - dY;
Y_H = Y_M + dY;
set(gca,'YLim',[Y_L,Y_H]);
set(gca,'YTick',[Y_L,Y_M,Y_H]);

if rmax == 5
% ### Subplot
subplot(rmax,cmax,2+cmax*3);
plot(t_vec,w3_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
dY = 0.003;
Y_M = 0.9994;
Y_L = Y_M - dY;
Y_H = Y_M + dY;
set(gca,'YLim',[Y_L,Y_H]);
set(gca,'YTick',[Y_L,Y_M,Y_H]);

% ### Subplot
subplot(rmax,cmax,2+cmax*4);
plot(t_vec,w8_vec,'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
dY = 0.003;
Y_M = 0.9994;
Y_L = Y_M - dY;
Y_H = Y_M + dY;
set(gca,'YLim',[Y_L,Y_H]);
set(gca,'YTick',[Y_L,Y_M,Y_H]);
end

xlabel({'Time (s)'},'interpreter','latex');

if enable_save; print(gcf,'Case_Composite_Simulation.png','-dpng','-r600'); end 
end