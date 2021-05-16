% This function plots the simulation results of the single-SG-infinite-bus
% system.

% Author(s): Yitong Li

%%
clear all
close all
clc

%%
enable_save = 0;

%%
% Load simulation results
load('Case_SG.mat');

% stop

% Organize data
fn = 0;
t_shift = -10;
t_vec = Case_SG.time;
t_vec = t_vec + t_shift;

% Set w to an array
vdq_vec     = Case_SG.signals(1).values;
idq_vec     = Case_SG.signals(2).values;
vabc_vec    = Case_SG.signals(3).values;
iabc_vec    = Case_SG.signals(4).values;
w_vec_      = Case_SG.signals(5).values;
theta_vec_  = Case_SG.signals(6).values;
pq_vec      = Case_SG.signals(7).values;
pqfil_vec   = Case_SG.signals(8).values;

N = length(w_vec_);
for i = 1:N
    w_vec(i) = w_vec_(1,1,i);
end

X_L = -1;
X_H = 10;
LineWidth = 1;
XTick = [-1,0,5,10];
XLim = [X_L,X_H];

%% Plot
if 1
    
fn = fn+1;
figure(fn)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.22 0.28]); % position of left-bottow cornor + length/depth of figure
plot(t_vec,w_vec,'LineWidth',LineWidth); grid on; hold on;
dYw = 0.002;
Yw_L = 1 - dYw;
Yw_H = 1 + dYw;
YLim = [Yw_L,Yw_H];
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',YLim);
xlabel({'Time (s)'},'interpreter','latex');
ylabel({'$\omega$ (pu)'},'interpreter','latex');

if enable_save
    print(gcf,'Case_SG_Simulation.png','-dpng','-r600'); 
end

end

%% Plot
if 0
fn = fn+1;
figure(fn)
% Set common
% set(gcf,'units','normalized','outerposition',[0.1 0.1 0.28 0.45]); % position of left-bottow cornor + length/depth of figure
subplot(1,2,1)
plot(t_vec,w_vec,'LineWidth',LineWidth); grid on; hold on;
dYw = 0.002;
Yw_L = 1 - dYw;
Yw_H = 1 + dYw;
set(gca,'XLim',[X_L,X_H]);
set(gca,'YLim',[Yw_L,Yw_H]);
xlabel({'Time (s)'},'interpreter','latex');
ylabel({'$\omega$ (pu)'},'interpreter','latex');

subplot(1,2,2)
plot(t_vec,-pqfil_vec(:,1),'LineWidth',LineWidth); grid on; hold on;
set(gca,'XLim',[X_L,X_H]);
% set(gca,'YLim',[Yw_L,Yw_H]);
xlabel({'Time (s)'},'interpreter','latex');
ylabel({'$p$ (pu)'},'interpreter','latex');

if enable_save
    print(gcf,'Case1_Simulation.png','-dpng','-r600'); 
end

end

