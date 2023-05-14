clear all
clc
close all

%%
RgbBlue = [0, 0.4470, 0.7410];
RgbRed = [0.8500, 0.3250, 0.0980];

PoleVq{1} = load('pole_ss_vq_150Hz').pole_ss;
PoleVq{2} = load('pole_ss_vq_210Hz').pole_ss;
PoleVq{3} = load('pole_ss_vq_300Hz').pole_ss;

PoleW{1} = load('pole_ss_W_150Hz').pole_ss;
PoleW{2} = load('pole_ss_W_210Hz').pole_ss;
PoleW{3} = load('pole_ss_W_300Hz').pole_ss;

%%
figure(1)
FigSize = [0.1 0.1 0.25 0.5];
set(gcf,'units','normalized','outerposition',FigSize);

for i = 1:length(PoleVq)
    PlotPole(PoleVq{i},'x',RgbBlue);
    PlotPole(PoleW{i},'o',RgbRed);
end
axis([-100,20,-150,150])
legend('vq-PLL','W-PLL')

print(gcf,['0 Sim_InfBusGfl\PolePLL.png'],'-dpng','-r600');

%%
function PlotPole(Pole,Type,Color)
    scatter(real(Pole),imag(Pole),Type,'LineWidth',1,'MarkerEdgeColor',Color); hold on; grid on;
    % scatter(real(Pole),imag(Pole),Type,'LineWidth',1); hold on; grid on;
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
end
