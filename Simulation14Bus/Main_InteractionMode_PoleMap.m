% Author(s): Yitong Li

clear all
clc
close all

Eig_Ref = load('Pole_Reference').EigVecHz;
Eig_DC30Hz = load('Pole_IBR8_DC30Hz').EigVecHz;
Eig_AC150Hz = load('Pole_IBR8_AC150Hz').EigVecHz;
Eig_PLL30Hz = load('Pole_IBR8_PLL30Hz').EigVecHz;

EnableSave = 0;

FigN = 0;
FigureSize = [0.1 0.1 0.17 0.38];

FigN = FigN+1;
figure(FigN);
set(gcf,'units','normalized','outerposition',FigureSize);
scatter(real(Eig_Ref),imag(Eig_Ref),'x','LineWidth',1.5); hold on; grid on;
scatter(real(Eig_DC30Hz),imag(Eig_DC30Hz),'x','LineWidth',1.5); hold on; grid on;
xlabel('Real Part (Hz)','interpreter','latex');
ylabel('Imaginary Part (Hz)','interpreter','latex')
axis([-40,10,-100,100]);
yticks([-100,-50,0,50,100]);
xticks([-40,-30,-20,-10,0,10]);
if EnableSave
    print(gcf,'Simulation/Figure/InteractionModePoleDC.png','-dpng','-r600');
end

FigN = FigN+1;
figure(FigN);
set(gcf,'units','normalized','outerposition',FigureSize);
scatter(real(Eig_Ref),imag(Eig_Ref),'x','LineWidth',1.5); hold on; grid on;
scatter(real(Eig_PLL30Hz),imag(Eig_PLL30Hz),'x','LineWidth',1.5); hold on; grid on;
xlabel('Real Part (Hz)','interpreter','latex');
ylabel('Imaginary Part (Hz)','interpreter','latex')
axis([-40,10,-100,100]);
yticks([-100,-50,0,50,100]);
xticks([-40,-30,-20,-10,0,10]);
if EnableSave
    print(gcf,'Simulation/Figure/InteractionModePolePLL.png','-dpng','-r600');
end

FigN = FigN+1;
figure(FigN);
set(gcf,'units','normalized','outerposition',FigureSize);
scatter(real(Eig_Ref),imag(Eig_Ref),'x','LineWidth',1.5); hold on; grid on;
scatter(real(Eig_AC150Hz),imag(Eig_AC150Hz),'x','LineWidth',1.5); hold on; grid on;
xlabel('Real Part (Hz)','interpreter','latex');
ylabel('Imaginary Part (Hz)','interpreter','latex')
axis([-40,10,-100,100]);
yticks([-100,-50,0,50,100]);
xticks([-40,-30,-20,-10,0,10]);
if EnableSave
    print(gcf,'Simulation/Figure/InteractionModePoleAC.png','-dpng','-r600');
end