clear all
clc
close all

load('Sim118BusResult.mat');

Time = out.App1{1}.Values.Time;
Time = Time - 4.3;

AllBus = [1:118];
FloatBus = [2,3,5,7,9,11,13,14,16,17,20:23,28:30,33,35,37:39,41,43:45,47,48,50:53,57,58,60,63,64,67,68,71,75,78,79,81:84,86,88,93:98,101,102,106,108,109,114,115,117,118];
AppIndex = setdiff(AllBus,FloatBus);

for m = 1:length(AppIndex)
    n = AppIndex(m);
    v_dq{n} = eval(['out.App' num2str(n) '{1}.Values.Data']);
    v_m{n} = ArrayNorm(v_dq{n});
    w_{n} = eval(['out.App' num2str(n) '{5}.Values.Data']);
    for k = 1:length(w_{n})
        w{n}(k) = w_{n}(:,:,k);
    end
    w{n} = transpose(w{n});
end

figure(1)
LineWidth = 1;
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.25 0.55]);

 for m = 1:length(AppIndex)
    
    n = AppIndex(m);
    
    subplot(2,1,1)
    plot(Time,w{n},'LineWidth',LineWidth); hold on; grid on;
    xlim([0,1.2]);
    ylim([0.96,1.02]);
    ylabel('Frequency (pu)')

    subplot(2,1,2)
    plot(Time,v_m{n},'LineWidth',LineWidth); hold on; grid on;
    xlim([0,1.2]);
    ylim([0.7,1.3]);
    ylabel('Voltage (pu)')
    xlabel('Time (s)')
 end
 
 if 1
 print(gcf,['Simulation118Bus\Sim118Bus_SimEmt.png'],'-dpng','-r600');
 end
 
%% 
function Output = ArrayNorm(Input)
    [r,c] = size(Input);
    for i = 1:r
        Output(i) = norm(Input(i,:));
    end
    Output = transpose(Output);
end