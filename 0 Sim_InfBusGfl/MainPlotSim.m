clear all
clc
close all

%%
out{1} = load('Emt_vqPLL10Hz_CC300Hz.mat').out;
out{2} = load('Emt_WPLL10Hz_CC300Hz.mat').out;

Time = out{1}.GflData{1}.Values.Time;
Timeshift = 0.16;
Time = Time - Timeshift;

for i = 1:length(out)
    vdq{i} = out{i}.GflData{1}.Values.Data;
    vm{i} = transpose(ArrayNorm(vdq{i}));
    idq{i} = out{i}.GflData{2}.Values.Data;
    im{i} = transpose(ArrayNorm(idq{i}));
    dtheta{i} = out{i}.GflData{6}.Values.Data(:,2);
    dtheta{i} = dtheta{i} - dtheta{i}(100);
end

% for i = 1:2
%     theta_ = theta{i};
%     for j = 1:length(theta{i})
%         theta__(j) = theta_(:,:,j);
%     end
%     theta{i} = theta__;
% end
% 
% for i = 1:2
%     dtheta{i} = theta{i} - theta{i}(100);
%     dtheta{i} = dtheta{i} + pi;
%     dtheta{i} = mod(dtheta{i},2*pi);
%     dtheta{i} = dtheta{i} - pi;
%     dtheta{i} = DisableMod(dtheta{i});
% end

%%
LineWidth1 = 1.6;
LineWidth2 = 1.2;
XLimit = [0,0.16];
XTicks = [0,0.04,0.08,0.12];

figure(1)
FigSize = [0.1 0.1 0.25 0.6];
set(gcf,'units','normalized','outerposition',FigSize);
% for i = 1:length(out)
%     subplot(3,1,1)
%     plot(Time,dtheta{i},'linewidth',LineWidth); grid on; hold on;
%     xlim(XLimit);
%     subplot(3,1,2)
%     plot(Time,vm{i},'linewidth',LineWidth); grid on; hold on;
%     xlim(XLimit);
%     subplot(3,1,3)
%     plot(Time,im{i},'linewidth',LineWidth); grid on; hold on;
%     xlim(XLimit);
% end

subplot(3,1,1)
plot(Time,dtheta{1},'linewidth',LineWidth1); grid on; hold on;
subplot(3,1,2)
plot(Time,vm{1},'linewidth',LineWidth1); grid on; hold on;
subplot(3,1,3)
plot(Time,im{1},'linewidth',LineWidth1); grid on; hold on;

subplot(3,1,1)
plot(Time,dtheta{2},'--','linewidth',LineWidth2); grid on; hold on;
xlim(XLimit);
xticks(XTicks);
ylim([-0.15,0.05])
yticks([-0.15,-0.1,-0.05,0,0.05])
ylabel('Angle (rad)')
subplot(3,1,2)
plot(Time,vm{2},'--','linewidth',LineWidth2); grid on; hold on;
xlim(XLimit);
xticks(XTicks);
ylim([0,1.5]);
ylabel('Voltage (pu)')
subplot(3,1,3)
plot(Time,im{2},'--','linewidth',LineWidth2); grid on; hold on;
xlim(XLimit);
xticks(XTicks);
ylim([0,3]);
ylabel('Current (pu)')
xlabel('Time (s)')

legend('vq-PLL','W-PLL')

print(gcf,['0 Sim_InfBusGfl\SimPLL.png'],'-dpng','-r600');

%%
function Output = ArrayNorm(Input)
    [r,c] = size(Input);
    for i = 1:r
        Output(i) = norm(Input(i,:));
    end
end

function thetaNew = DisableMod(theta)
    len = length(theta);
    theta(len+1) = theta(len);
    for i = 1:len
        if theta(i+1) - theta(i) > pi
            theta(i+1:end) = theta(i+1:end) - 2*pi;
        elseif theta(i+1) - theta(i) < -pi
            theta(i+1:end) = theta(i+1:end) + 2*pi;
        end
    end
    thetaNew = theta(1:len);
end
