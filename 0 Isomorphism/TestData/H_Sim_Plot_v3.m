% This script plots the inertia and bandwidth analysis

clear all
clc
close all

mfile_name = mfilename('fullpath');
[RootPath,~,~]  = fileparts(mfile_name);
cd(RootPath);

%% Enables
Enable_SaveFigure = 1;

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
p_SG = -pq_SG(:,1);
s_SG = pq_SG(:,1) + 1i*pq_SG(:,2);
s_SG_abs = abs(s_SG);
s_SG_int = zeros(1,TimeIndex-1);
p_SG_int = zeros(1,TimeIndex-1);
for i = TimeIndex:TimeLength
    s_SG_int(i) = s_SG_int(i-1)+(s_SG_abs(i) - s_SG_abs(TimeLength))*DeltaTime;
    p_SG_int(i) = p_SG_int(i-1)+(p_SG(i) - p_SG(TimeLength))*DeltaTime;
end

vdq_IBR = H_SimData_IBR.signals(1).values;
w_IBR_ = H_SimData_IBR.signals(5).values;
for i = 1:length(w_IBR_)
w_IBR(i) = w_IBR_(:,:,i);
end
pq_IBR = H_SimData_IBR.signals(7).values;
p_IBR = -pq_IBR(:,1);
s_IBR = pq_IBR(:,1) + 1i*pq_IBR(:,2);
s_IBR_abs = abs(s_IBR);
s_IBR_int = zeros(1,TimeIndex-1);
p_IBR_int = zeros(1,TimeIndex-1);
for i = TimeIndex:TimeLength
    s_IBR_int(i) = s_IBR_int(i-1)+(s_IBR_abs(i) - s_IBR_abs(TimeLength))*DeltaTime;
    p_IBR_int(i) = p_IBR_int(i-1)+(p_IBR(i) - p_IBR(TimeLength))*DeltaTime;
end

%% Plot Settings
FigNum = 0;
LineWidth = 1;

FigRowMax = 2;
FigColumnMax = 1;

Time = Time - 1.8;
t_Limit = [0,0.6];
t_Ticks = [0,0.2,0.4,0.6];
w_Limit = [0.98,1.02]; 
w_Ticks = [0.98,1,1.02];
s_Limit = [0,1.5];
s_Ticks = [0,0.5,1,1.5];
% e_Limit = [-15,15]*1e-3;
% e_Ticks = [-15,0,15]*1e-3;
e_Limit = [-0.2,1];
e_Ticks = [-0.2,0,0.2,0.4,0.6,0.8,1];

FigSize = [0.1 0.1 0.18 0.5];

ColorLower = [0,1,1]; ColorUpper = [0,0,1];       % light blue to dark blue

F0 =60;
FreqLower = 0.4*F0;
FreqUpper = 2*F0;
fb(1) = 25.239723972397233;
fb(2) = 98.411041104110396;
for k = 1:length(fb)
FreqFactor = (fb(k)-FreqLower)/(FreqUpper-FreqLower);
LineColor{k} = (ColorUpper - ColorLower) * FreqFactor + ColorLower;
end

%% Plot
% SG
FigNum = FigNum + 1;
figure(FigNum)
set(gcf,'units','normalized','outerposition',FigSize);

subplot(FigRowMax,FigColumnMax,1)
plot(Time,s_SG_abs,'LineWidth',LineWidth,'Color',LineColor{1}); grid on; hold on;
ylabel('$|S|$ (pu)','interpreter','latex')
xlim(t_Limit);
xticks(t_Ticks)
ylim(s_Limit);
yticks(s_Ticks);

subplot(FigRowMax,FigColumnMax,2)
plot(Time,abs(p_SG_int)*60,'LineWidth',LineWidth,'Color',LineColor{1}); grid on; hold on;
% xlabel('Time (s)','interpreter','latex')
xlabel('Time (s)')
ylabel('$E$ (pu)','interpreter','latex')
xlim(t_Limit);
xticks(t_Ticks)
ylim(e_Limit);
yticks(e_Ticks);

if Enable_SaveFigure
    print(gcf,'H_Sim_SG.png','-dpng','-r600');
end

% IBR
FigNum = FigNum + 1;
figure(FigNum)
set(gcf,'units','normalized','outerposition',FigSize);
subplot(FigRowMax,FigColumnMax,1)

subplot(FigRowMax,FigColumnMax,1)
plot(Time,s_IBR_abs,'LineWidth',LineWidth,'Color',LineColor{2}); grid on; hold on;
ylabel('$|S|$ (pu)','interpreter','latex')
xlim(t_Limit);
xticks(t_Ticks)
ylim(s_Limit);
yticks(s_Ticks);

subplot(FigRowMax,FigColumnMax,2)
plot(Time,abs(p_IBR_int)*60,'LineWidth',LineWidth,'Color',LineColor{2}); grid on; hold on;
% xlabel('Time (s)','interpreter','latex')
xlabel('Time (s)')
ylabel('$E$ (pu)','interpreter','latex')
xlim(t_Limit);
xticks(t_Ticks)
ylim(e_Limit);
yticks(e_Ticks);

if Enable_SaveFigure
    print(gcf,'H_Sim_IBR.png','-dpng','-r600');
end

% IBR zoom-in plot
FigNum = FigNum + 1;
figure(FigNum)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.11 0.32]);
subplot(FigRowMax,FigColumnMax,1)

subplot(FigRowMax,FigColumnMax,1)
plot(Time,s_IBR_abs,'LineWidth',LineWidth,'Color',LineColor{2}); grid on; hold on;
xlim([0.1,0.4]);
xticks([0.1,0.2,0.3,0.4])
ylim([0.45,0.55]);
yticks([0.45,0.5,0.55]);

subplot(FigRowMax,FigColumnMax,2)
plot(Time,p_IBR_int*60,'LineWidth',LineWidth,'Color',LineColor{2}); grid on; hold on;
xlim([0.1,0.4]);
xticks([0.1,0.2,0.3,0.4])
ylim([-1,1]*1e-2);
yticks([-1,0,1]*1e-2);

if Enable_SaveFigure
    print(gcf,'H_Sim_IBR_ZoomIn.png','-dpng','-r600');
end