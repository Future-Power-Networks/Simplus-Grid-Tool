% Author(s): Yitong Li

clear all
clc
close all

EnableSave = 1;

FigN = 0;

load('SimSwingMode');
SgOmega = out.SgOmega;

time = SgOmega{1}.Values.Time;
w1_ = SgOmega{1}.Values.Data;
w2_ = SgOmega{2}.Values.Data;
w3_ = SgOmega{3}.Values.Data;

for i = 1:length(w1_)
    w1(i) = w1_(:,:,i);
end
for i = 1:length(w2_)
    w2(i) = w2_(:,:,i);
end
for i = 1:length(w1_)
    w3(i) = w3_(:,:,i);
end

FigN = FigN+1;
figure(FigN)

time = time - 5.5;
LineWidth = 1;
XLim = [0,3];
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.25 0.45]);

subplot(3,1,1)
plot(time,w1,'LineWidth',LineWidth); grid on; hold on;
ylabel('$\omega_1$ (pu)','interpreter','latex')
xlim(XLim);
ylim([0.994,1.006]);
% xticks(XTicks);
% ylim([0.999,1.001]);
% yticks([0.999,1,1.001]);

subplot(3,1,2)
plot(time,w2,'LineWidth',LineWidth); grid on; hold on;
ylabel('$\omega_3$ (pu)','interpreter','latex')
xlim(XLim);
ylim([0.994,1.006]);

subplot(3,1,3)
plot(time,w3,'LineWidth',LineWidth); grid on; hold on;
ylabel('$\omega_6$ (pu)','interpreter','latex')
xlim(XLim);
ylim([0.99,1.01]);

xlabel('Time (s)','interpreter','latex')

if EnableSave
    print(gcf,'Simulation/Figure/SwingModeEmt.png','-dpng','-r600');
end