% Author(s): Yitong Li

clear all
clc
close all

%%
Eig_Ref = load('EigVecHz_118Bus_Ref.mat').EigVecHz;
Eig_Bv_Large = load('EigVecHz_118Bus_Bv_1e-1.mat').EigVecHz;
Eig_Gv_Large = load('EigVecHz_118Bus_Gv_1e-1.mat').EigVecHz;

%%
figure(1)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.28 0.55]);
PlotPoleMap(Eig_Ref);
PlotPoleMap(Eig_Bv_Large); 
PlotPoleMap(Eig_Gv_Large);
axis([-60,10,-150,150]);
title('Eigenvalue')

legend('Descriptor or Cv=1e-5 or Rv=1/1e-5','Cv=1e-1','Rv=1/1e-1','location','southeast')

if 1
print(gcf,['Simulation118Bus\118Bus_PoleMap.png'],'-dpng','-r600');
end

figure(2)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.18 0.55]);
PlotPoleMap(Eig_Ref);
PlotPoleMap(Eig_Bv_Large);
PlotPoleMap(Eig_Gv_Large);

axis([-2,0,28,28.6]);
title('Zoom In')

if 1
print(gcf,['Simulation118Bus\118Bus_PoleMap_ZoomIn.png'],'-dpng','-r600');
end

%%
function PlotPoleMap(Eig)
    
    scatter(real(Eig),imag(Eig),'x','LineWidth',1.5); hold on; grid on;
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
    
end