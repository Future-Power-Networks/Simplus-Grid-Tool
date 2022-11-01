clear all
close all
clc

cd(fileparts(mfilename('fullpath')));

% GFL analysis
if 0
for i = 2:2
UserDataType = 1;
UserDataName = ['GflTest_v',num2str(i)];
SimplusGT.Toolbox.Main();
% save([UserDataName,'.mat'],'EigVecHz')
PlotPoleMap(EigVecHz,9999)
clear all
end
end

% GFM analysis
if 1
for i = 1:1
UserDataType = 1;
UserDataName = ['GfmTest_v',num2str(i)];
SimplusGT.Toolbox.Main();
% save([UserDataName,'.mat'],'EigVecHz')
PlotPoleMap(EigVecHz,9999)
clear all
end
end

%%
function PlotPoleMap(EigVecHz,FigN)
  
    figure(FigN);
    
    subplot(1,2,1)
    scatter(real(EigVecHz),imag(EigVecHz),'x','LineWidth',1.5); hold on; grid on;
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
    title('Global pole map');
    
	subplot(1,2,2)
    scatter(real(EigVecHz),imag(EigVecHz),'x','LineWidth',1.5); hold on; grid on;
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
    title('Zoomed pole map');
    axis([-80,20,-150,150]);
end