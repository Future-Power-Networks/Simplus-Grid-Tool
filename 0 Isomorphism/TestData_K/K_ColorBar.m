clear all
clc
close all

ColorStepSize = 100;
ColorLower = [1,1,1];
ColorUpper = [1,0.5,0.5];
GradRed     = linspace(ColorLower(1),ColorUpper(1),ColorStepSize)';
GradGreen   = linspace(ColorLower(2),ColorUpper(2),ColorStepSize)';
GradBlue    = linspace(ColorLower(3),ColorUpper(3),ColorStepSize)';
colormap([GradRed GradGreen GradBlue]);
colorbar();
caxis([0 0.25]);

print(gcf,['Graph_Impedance_ColorBar.png'],'-dpng','-r600');
