clear all
clc
close all

load('Sim68BusResult.mat');

Time = out.App1{1}.Values.Time;
Time = Time - 4.3;

for i = 1:16
    v_dq{i} = eval(['out.App' num2str(i) '{1}.Values.Data']);
    v_m{i} = ArrayNorm(v_dq{i});
    w_{i} = eval(['out.App' num2str(i) '{5}.Values.Data']);
    for k = 1:length(w_{i})
        w{i}(k) = w_{i}(:,:,k);
    end
    w{i} = transpose(w{i});
end

figure(1)
LineWidth = 1;
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.25 0.55]);

 for i = 1:16
    subplot(2,1,1)
    plot(Time,w{i},'LineWidth',LineWidth); hold on; grid on;
    xlim([0,1.2]);
    ylim([0.96,1.02]);
    ylabel('Frequency (pu)')

    subplot(2,1,2)
    plot(Time,v_m{i},'LineWidth',LineWidth); hold on; grid on;
    xlim([0,1.2]);
    ylim([0.7,1.3]);
    ylabel('Voltage (pu)')
    xlabel('Time (s)')
 end
 
 print(gcf,['Simulation68Bus\Sim68Bus_SimEmt.png'],'-dpng','-r600');

%% 
function Output = ArrayNorm(Input)
    [r,c] = size(Input);
    for i = 1:r
        Output(i) = norm(Input(i,:));
    end
    Output = transpose(Output);
end