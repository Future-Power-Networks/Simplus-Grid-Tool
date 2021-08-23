% This script plots the H analysis

clear all
clc
close all


%% Load data
H_SimData_SG = load('H_SimData_SG').H_SimData_SG;
H_SimData_IBR = load('H_SimData_IBR').H_SimData_IBR;

Time = H_SimData_SG.time;
DeltaTime = Time(2)-Time(1);
TimeIndex = find(Time == 2);
TimeLength = length(Time);

vdq_SG = H_SimData_SG.signals(1).values;
w_SG_ = H_SimData_SG.signals(5).values;
for i = 1:length(w_SG_)
w_SG(i) = w_SG_(:,:,i);
end
pq_SG = H_SimData_SG.signals(7).values;
s_SG = pq_SG(:,1) + 1i*pq_SG(:,2);
s_SG_abs = abs(s_SG);
s_SG_int = zeros(1,TimeIndex-1);
for i = TimeIndex:TimeLength
    s_SG_int(i) = s_SG_int(i-1)+(s_SG_abs(i) - s_SG_abs(TimeLength))*DeltaTime;
end

vdq_IBR = H_SimData_IBR.signals(1).values;
w_IBR_ = H_SimData_IBR.signals(5).values;
for i = 1:length(w_IBR_)
w_IBR(i) = w_IBR_(:,:,i);
end
pq_IBR = H_SimData_IBR.signals(7).values;
s_IBR = pq_IBR(:,1) + 1i*pq_IBR(:,2);
s_IBR_abs = abs(s_IBR);
s_IBR_int = zeros(1,TimeIndex-1);
for i = TimeIndex:TimeLength
    s_IBR_int(i) = s_IBR_int(i-1)+(s_IBR_abs(i) - s_IBR_abs(TimeLength))*DeltaTime;
end

%% Plot Settings
FigNum = 0;
LineWidth = 1.2;

FigRowMax = 3;
FigColumnMax = 2;

Time = Time - 1.8;
t_Limit = [0,0.8];
t_Ticks = [0,0.4,0.8];
w_Limit = [0.98,1.01]; 
s_Limit = [0,1];
s_Ticks = [0,0.5,1];
e_Limit = [-5,10]*1e-3;
e_Ticks = [-5,0,5,10]*1e-3;

%% Plot
FigNum = FigNum + 1;
figure(FigNum)
subplot(FigRowMax,FigColumnMax,1)
plot(Time,w_SG,'LineWidth',LineWidth); grid on; hold on;
%xlabel('Time (s)','interpreter','latex')
ylabel('$\omega$ (pu)','interpreter','latex')
xlim(t_Limit);
xticks(t_Ticks)
ylim(w_Limit);
% yticks(w_Ticks);

subplot(FigRowMax,FigColumnMax,2)
plot(Time,w_IBR,'LineWidth',LineWidth); grid on; hold on;
%xlabel('Time (s)','interpreter','latex')
ylabel('$\omega$ (pu)','interpreter','latex')
xlim(t_Limit);
xticks(t_Ticks)
ylim(w_Limit);

subplot(FigRowMax,FigColumnMax,3)
plot(Time,s_SG_abs,'LineWidth',LineWidth); grid on; hold on;
xlabel('Time (s)','interpreter','latex')
ylabel('$S$ (pu)','interpreter','latex')
xlim(t_Limit);
xticks(t_Ticks)
ylim(s_Limit);
yticks(s_Ticks);

subplot(FigRowMax,FigColumnMax,4)
plot(Time,s_IBR_abs,'LineWidth',LineWidth); grid on; hold on;
xlabel('Time (s)','interpreter','latex')
ylabel('$S$ (pu)','interpreter','latex')
xlim(t_Limit);
xticks(t_Ticks)
ylim(s_Limit);
yticks(s_Ticks);

subplot(FigRowMax,FigColumnMax,5)
plot(Time,s_SG_int,'LineWidth',LineWidth); grid on; hold on;
xlabel('Time (s)','interpreter','latex')
ylabel('$E$ (pu)','interpreter','latex')
xlim(t_Limit);
xticks(t_Ticks)
ylim(e_Limit);
yticks(e_Ticks);

subplot(FigRowMax,FigColumnMax,6)
plot(Time,s_IBR_int,'LineWidth',LineWidth); grid on; hold on;
xlabel('Time (s)','interpreter','latex')
ylabel('$E$ (pu)','interpreter','latex')
xlim(t_Limit);
xticks(t_Ticks)
ylim(e_Limit);
yticks(e_Ticks);