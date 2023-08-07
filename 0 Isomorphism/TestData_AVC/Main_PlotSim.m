% Plot the simulation results
clear all
clc
close all

DataName = 'out_AVC';
% DataName = 'out_Iref';

TimeShift = 2;
Enable_SaveFigure = 1;
FigSize = [0.1 0.1 0.25 0.6];

load(DataName);

Time = out.Data_App1{6}.Values.Time;
Time = Time - TimeShift;

Fbus = [19,22,30,31,32,34,35,37,38,43,54,57,58,62,63,65,66];

NumRef = 13;
NumBus = 16;


for i = 1:NumBus
    if isempty(find(Fbus == i,1))
        theta{i} = eval(['out.Data_App' num2str(i) '{6}.Values.Data']);
        vdq{i} = eval(['out.Data_App' num2str(i) '{1}.Values.Data']);
        vm{i} = ArrayNorm(vdq{i});
        idq{i} = eval(['out.Data_App' num2str(i) '{2}.Values.Data']);
        im{i} = ArrayNorm(idq{i});
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
set(gcf,'units','normalized','outerposition',FigSize);
for i = 1:NumBus
    if isempty(find(Fbus == i,1))
        subplot(3,1,1)
        plot(Time,dtheta{i}); hold on; grid on;
        xlim([0,1.5]);
        % ylim([0,1]); 
        ylim([-1,2]); 
        ylabel('Angle (rad)');
        
        subplot(3,1,2)
        plot(Time,vm{i}); hold on; grid on;
        xlim([0,1.5]);
      	% ylim([0,1.5]);
        ylim([0,2]);
        ylabel('Voltage (pu)')
        
        subplot(3,1,3)
        plot(Time,im{i}); hold on; grid on;
        xlim([0,1.5]);
      	ylim([0,40]);
        ylabel('Current (pu)')
        
        xlabel('Time (s)')
    end
end

if Enable_SaveFigure
    print(gcf,['TestData_AVC/K_Sim_' DataName '.png'],'-dpng','-r600');
end