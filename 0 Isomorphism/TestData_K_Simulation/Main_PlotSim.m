% Plot the simulation results
clear all
clc
close all

% DataName = '0SG_Load';
% DataName = '0SG';
% DataName = '2SG';
DataName = '3SG';

TimeShift = 2;
Enable_SaveFigure = 1;
FigSize = [0.1 0.1 0.35 0.5];

load(['out_' DataName]);

Time = out.Data_App1{6}.Values.Time;
Time = Time - TimeShift;

Fbus = [19,22,30,31,32,34,35,37,38,43,54,57,58,62,63,65,66];


NumRef = 17;
NumBus = 68;

if strcmp(DataName,'0SG_Load')
    NumRef = 13;
    NumBus = 16;
end


for i = 1:NumBus
    if isempty(find(Fbus == i,1))
        theta{i} = eval(['out.Data_App' num2str(i) '{6}.Values.Data']);
        vdq{i} = eval(['out.Data_App' num2str(i) '{1}.Values.Data']);
        vm{i} = ArrayNorm(vdq{i});
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
        subplot(2,1,1)
        plot(Time,dtheta{i}); hold on; grid on;
        xlim([0,1.5]);
        if strcmp(DataName,'0SG_Load') || strcmp(DataName,'0SG') || strcmp(DataName,'3SG')
            ylim([0,1]);        % 0SG, 3SG, 0SG_Load
        elseif strcmp(DataName,'2SG')
            ylim([-5,30]);     % 2SG
        else
            error(['Error!']);
        end
        ylabel('Angle (rad)');
        subplot(2,1,2)
        plot(Time,vm{i}); hold on; grid on;
        xlim([0,1.5]);
        if strcmp(DataName,'0SG')
            ylim([0,4]);      % 0SG
        elseif strcmp(DataName,'2SG')
            ylim([0,8]);        % 2SG
        elseif strcmp(DataName,'3SG') || strcmp(DataName,'0SG_Load')
            ylim([0,1.5]);        % 3SG, 0SG_Load
        else
            error(['Error!']);
        end
        ylabel('Voltage (pu)')
        xlabel('Time (s)')
    end
end

if Enable_SaveFigure
    print(gcf,['K_Sim_' DataName '.png'],'-dpng','-r600');
end