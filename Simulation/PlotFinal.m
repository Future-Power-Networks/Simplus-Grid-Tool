% Plot the simulation results
%
% Author(s): Yitong Li

clear all
clc
close all

%% Select data
% DataName = 'LocalMode_HighInertia';
% DataName = 'LocalMode_MedInertia';
% DataName = 'LocalMode_LowInertia';

DataName = 'InterAreaMode_HighInertia';
% DataName = 'InterAreaMode_MedInertia';
% DataName = 'InterAreaMode_LowInertia';



%% Plot
TimeShift = 2;
if strcmp(DataName,'InterAreaMode_LowInertia') || strcmp(DataName,'LocalMode_LowInertia') || ...
   strcmp(DataName,'InterAreaMode_MedInertia') || strcmp(DataName,'LocalMode_MedInertia')
    TimeShift = 8;
end
Enable_SaveFigure = 1;
% FigSize = [0.1 0.1 0.35 0.7];
FigSize = [0.1 0.1 0.25 0.6];

load(['out_' DataName]);

Time = out.Data_App1{6}.Values.Time;
Time = Time - TimeShift;

Fbus = [19,22,30,31,32,34,35,37,38,43,54,57,58,62,63,65,66];

NumRef = 13;
NumBus = 16;

for i = 1:NumBus
    if isempty(find(Fbus == i,1))
        w_{i} = eval(['out.Data_App' num2str(i) '{5}.Values.Data']);
        theta{i} = eval(['out.Data_App' num2str(i) '{6}.Values.Data']);
        vdq{i} = eval(['out.Data_App' num2str(i) '{1}.Values.Data']);
        vm{i} = ArrayNorm(vdq{i});
    end
end
for i = 1:NumBus
    for k = 1:length(w_{i})
        w{i}(k) = w_{i}(:,:,k);
    end
end

for i = 1:NumBus
    if isempty(find(Fbus == i,1))
        len = length(theta{i});
        theta_ = theta{i};
        for j = 1:len
         	theta__(j) = theta_(:,:,j);
        end
        theta{i} = theta__;
    end
end

% IndexMinTheta = 1;
% MinTheta = theta{1}(1);
% for i = 2:NumBus
%     theta_a = mod(theta{i-1}(1)+pi, 2*pi)-pi;
%     theta_b = mod(theta{i}(1)+pi, 2*pi)-pi;
%     if theta_a > theta_b
%         MinTheta = theta{i}(1);
%         IndexMinTheta = i;
%     end
% end

for i = 1:NumBus
    if isempty(find(Fbus == i,1))
        dtheta{i} = theta{i} - theta{NumRef};
        dtheta{i} = dtheta{i} + pi;
        dtheta{i} = mod(dtheta{i},2*pi);
        dtheta{i} = dtheta{i} - pi;
        dtheta{i} = DisableMod(dtheta{i});
    end
end

figure(1)
LineWidth = 1;
set(gcf,'units','normalized','outerposition',FigSize);
for i = 1:NumBus
% for i = 1
    if isempty(find(Fbus == i,1))
        subplot(3,1,1)
        plot(Time,w{i}(:,:,1),'LineWidth',LineWidth); hold on; grid on;
        xlim([0,2.5]);
        % ylim([0.88,1.12]);      % Local mode: low, medieum
        % ylim([0.95,1.1]);     	% Inter-area mode: low, medieum
        % ylim([0.95,1.2]);           % Local mode: high
        %                           % Inter-area mode: high
        ylabel('Frequency (pu)')
        subplot(3,1,2)
        plot(Time,dtheta{i}*180/pi,'LineWidth',LineWidth); hold on; grid on;
        % ylim([-20,400]);        % Local mode: high
        %                       % Inter mode: high
        % ylim([-10,120]);        % All modes: low, medieum
        xlim([0,2.5]);
        ylabel('Angle (Degree)');
        subplot(3,1,3)
        plot(Time,vm{i},'LineWidth',LineWidth); hold on; grid on;
        xlim([0,2.5]);
        ylim([0,1.5]);
        ylabel('Voltage (pu)')
        xlabel('Time (s)')
    end
end

if Enable_SaveFigure
    print(gcf,[DataName '.png'],'-dpng','-r600');
end