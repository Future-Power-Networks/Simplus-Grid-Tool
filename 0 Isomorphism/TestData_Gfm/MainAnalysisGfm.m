clear all
clc
close all

%%
RgbBlue = [0, 0.4470, 0.7410];
RgbRed = [0.8500, 0.3250, 0.0980];
RgbYellow = [0.9290, 0.6940, 0.1250];

%%
Pole_FullOrder = load('pole_sys_FullOrder').pole_sys;
Pole_NoShift_NoL = load('pole_T1cl_NoL').pole_T1cl;
Pole_NoShift_WithL = load('pole_T1cl_WithL').pole_T1cl;
Pole_WithShift_NoL = load('pole_T12cl_NoL').pole_T12cl;
Pole_WithShift_WithL = load('pole_T12cl_WithL').pole_T12cl;

%%
LineWdith = 1;
FigSize = [0.1,0.1,0.35,0.6];

%%
figure(1)
set(gcf,'units','normalized','outerposition',FigSize);
scatter(real(Pole_NoShift_WithL),imag(Pole_NoShift_WithL),'x','LineWidth',LineWdith,'MarkerEdgeColor',RgbYellow); hold on; grid on;
scatter(real(Pole_WithShift_WithL),imag(Pole_WithShift_WithL),'x','LineWidth',LineWdith,'MarkerEdgeColor',RgbBlue); hold on; grid on;
scatter(real(Pole_FullOrder),imag(Pole_FullOrder),'o','LineWidth',LineWdith,'MarkerEdgeColor',RgbRed); hold on; grid on;

axis([-0.015,0.005,-3,3]);
legend('No Channel Dynamics','Isomorphism','Full Order','location','southeast')
ylabel('Imaginary Part (Hz)')
xlabel('Real Part (Hz)')

print(gcf,['TestData_Gfm\PoleGfmFreqShift.png'],'-dpng','-r600');