% Author(s): Yitong Li

clc
close all

Eig_Ref = load('Pole_Reference').EigVecHz;

Eig_Sg6_D1 = load('Pole_Sg6_D1').EigVecHz;
Eig_Sg6_D20 = load('Pole_Sg6_D20').EigVecHz;

Eig_Sg13_D1 = load('Pole_Sg1Sg3_D1').EigVecHz;
Eig_Sg13_D20 = load('Pole_Sg1Sg3_D20').EigVecHz;

EnableSave = 1;

FigN = 0;

FigN = FigN+1;
figure(FigN);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.17 0.38]);
scatter(real(Eig_Sg6_D20),imag(Eig_Sg6_D20),'x','LineWidth',1.5); hold on; grid on;
scatter(real(Eig_Ref),imag(Eig_Ref),'x','LineWidth',1.5); hold on; grid on;
scatter(real(Eig_Sg6_D1),imag(Eig_Sg6_D1),'x','LineWidth',1.5); hold on; grid on;
axis([-1,0.5,-6,6]);
xticks([-1,-0.5,0,0.5]);
yticks([-6,-3,0,3,6]);
xlabel('Real Part (Hz)','interpreter','latex');
ylabel('Imaginary Part (Hz)','interpreter','latex')
if EnableSave
    print(gcf,'Simulation/Figure/SwingModePoleSg6.png','-dpng','-r600');
end

FigN = FigN+1;
figure(FigN);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.17 0.38]);
scatter(real(Eig_Sg13_D20),imag(Eig_Sg13_D20),'x','LineWidth',1.5); hold on; grid on;
scatter(real(Eig_Ref),imag(Eig_Ref),'x','LineWidth',1.5); hold on; grid on;
scatter(real(Eig_Sg13_D1),imag(Eig_Sg13_D1),'x','LineWidth',1.5); hold on; grid on;
axis([-1,0.5,-6,6]);
xticks([-1,-0.5,0,0.5]);
yticks([-6,-3,0,3,6]);
xlabel('Real Part (Hz)','interpreter','latex');
ylabel('Imaginary Part (Hz)','interpreter','latex')
if EnableSave
    print(gcf,'Simulation/Figure/SwingModePoleSg1Sg3.png','-dpng','-r600');
end
