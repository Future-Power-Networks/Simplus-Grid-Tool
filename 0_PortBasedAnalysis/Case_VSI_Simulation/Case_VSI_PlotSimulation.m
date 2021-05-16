% This function plots the simulation results of the single-SG-infinite-bus
% system.

% Author(s): Yitong Li

%%
clear all
close all

%%
enable_save = 0;

%%
% Load simulation results
load('Case_VSI.mat');

% Organize data
fn = 0;
t_shift = -1;
t_vec = Case_VSI.time;
t_vec = t_vec + t_shift;

% Set w to an array
vdq_vec     = Case_VSI.signals(1).values;
idq_vec     = Case_VSI.signals(2).values;
vabc_vec    = Case_VSI.signals(3).values;
iabc_vec    = Case_VSI.signals(4).values;
w_vec_      = Case_VSI.signals(5).values;
theta_vec_  = Case_VSI.signals(6).values;
pq_vec      = Case_VSI.signals(7).values;

N = length(w_vec_);
for i = 1:N
    w_vec(i) = w_vec_(1,1,i);
end

X_L = -0.05;
X_H = 0.2;
LineWidth = 1;
XTick = [-0.05,0,0.1,0.2];
XLim = [X_L,X_H];

%% Plot
if 1
    
fn = fn+1;
figure(fn)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.22 0.28]); % position of left-bottow cornor + length/depth of figure
plot(t_vec,w_vec,'LineWidth',LineWidth); grid on; hold on;
dYw = 0.2;
Yw_L = 1 - dYw;
Yw_H = 1 + dYw;
YLim = [Yw_L,Yw_H];
set(gca,'XLim',XLim);
set(gca,'XTick',XTick);
set(gca,'YLim',YLim);
xlabel({'Time (s)'},'interpreter','latex');
ylabel({'$\omega$ (pu)'},'interpreter','latex');

if enable_save
    print(gcf,'Case_VSI_Simulation.png','-dpng','-r600'); 
end

end