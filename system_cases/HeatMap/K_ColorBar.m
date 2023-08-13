% ColorStepSize = 100;
% ColorLower = [0,0,1];
% ColorUpper = [1,0,0];
% GradRed     = logspace(ColorLower(1),ColorUpper(1),ColorStepSize)'/10;
% GradGreen   = logspace(ColorLower(2),ColorUpper(2),ColorStepSize)'/10;
% GradBlue    = logspace(ColorLower(3),ColorUpper(3),ColorStepSize)'/10;
% colormap([ GradRed GradGreen GradBlue]);
% %colormap(flipud(jet))
% colorbar();
% caxis([0 100]);
% 
% print(gcf,['Graph_Impedance_ColorBar.png'],'-dpng','-r600');


c1 = [168, 0, 0]/255;
c2 = [214, 102, 100]/255;
c3 = [245, 216, 216]/255;
c4 = [255, 255, 255]/255;

for i=1:3
    g1=linspace(c1(i),c2(i),10);
    g2=linspace(c2(i),c3(i),15);
    g3=linspace(c3(i),c4(i),75);
    gx(:,i)=[g1,g2,g3]';
end

colormap(gx);
%colormap(flipud(jet))
colorbar();
caxis([0 10]);

print(gcf,['Graph_Impedance_ColorBar.png'],'-dpng','-r600');
