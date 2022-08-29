% Author(s): Yitong Li

clc
close all

EnableSave = 1;

FigN = 0;

FigN = FigN+1;
figure(FigN);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.35 0.45]);

subplot(1,2,1)
scatter(real(EigVecHz),imag(EigVecHz),'x','LineWidth',1.5); hold on; grid on;
xlabel('Real Part (Hz)','interpreter','latex');
ylabel('Imaginary Part (Hz)','interpreter','latex');
axis([-40,10,-100,100]);
yticks([-100,-50,0,50,100]);
xticks([-40,-30,-20,-10,0,10]);

subplot(1,2,2)
scatter(real(EigVecHz),imag(EigVecHz),'x','LineWidth',1.5); hold on; grid on;
axis([-1,0.5,-6,6]);
yticks([-6,-3,0,3,6]);
xticks([-1,-0.5,0,0.5]);

if EnableSave
    print(gcf,'Simulation/Figure/PoleMapRaw.png','-dpng','-r600');
end    

    