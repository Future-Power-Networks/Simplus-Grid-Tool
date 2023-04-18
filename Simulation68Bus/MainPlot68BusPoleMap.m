clear all
clc
close all

%%
Eig_Ref = load('EigVecHz_Ref.mat').EigVecHz;

Eig_Cv_Normal = load('EigVecHz_Cv_1e-3.mat').EigVecHz;
Eig_Cv_Small = load('EigVecHz_Cv_1e-12.mat').EigVecHz;
Eig_Cv_Large = load('EigVecHz_Cv_1e-1.mat').EigVecHz;

Eig_Gv_Normal = load('EigVecHz_Gv_1e-3.mat').EigVecHz;
Eig_Gv_Small = load('EigVecHz_Gv_1e-12.mat').EigVecHz;
Eig_Gv_Large = load('EigVecHz_Gv_1e-1.mat').EigVecHz;

%%
EnableSave = 1;

%%
figure(1)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.28 0.55]);
PlotPoleMap(Eig_Ref);
% PlotPoleMap(Eig_Cv_Normal);
% PlotPoleMap(Eig_Gv_Normal);

PlotPoleMap(Eig_Cv_Small); 
PlotPoleMap(Eig_Cv_Large); 
PlotPoleMap(Eig_Gv_Small);
PlotPoleMap(Eig_Gv_Large);

axis([-60,10,-150,150]);
title('Eigenvalue')

legend('Descriptor or Cv=1e-3 or Rv=1/1e-3','Cv=1e-12','Cv=1e-1','Rv=1/1e-12','Rv=1/1e-1','location','southeast')

if EnableSave == 1
print(gcf,['Simulation68Bus\68Bus_PoleMap.png'],'-dpng','-r600');
end

figure(2)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.18 0.55]);
PlotPoleMap(Eig_Ref);
% PlotPoleMap(Eig_Cv_Normal);
% PlotPoleMap(Eig_Gv_Normal);

PlotPoleMap(Eig_Cv_Small);
PlotPoleMap(Eig_Cv_Large);
PlotPoleMap(Eig_Gv_Small);
PlotPoleMap(Eig_Gv_Large);

axis([-5,5,25,58]);
title('Zoom In')

if EnableSave == 1
print(gcf,['Simulation68Bus\68Bus_PoleMap_ZoomIn.png'],'-dpng','-r600');
end

%%
function PlotPoleMap(Eig)
    
    scatter(real(Eig),imag(Eig),'x','LineWidth',1.5); hold on; grid on;
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
    
end