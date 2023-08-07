% This function is used to plot the wb-lambda contour: i.e., the
% trace of lambda when fixing wb to different values.

% Author(s): Yitong Li, Yunjie Gu

%%
clear all
clc
close all

mfile_name = mfilename('fullpath');
[RootPath,~,~]  = fileparts(mfile_name);
cd(RootPath);

ColorRGB();

%% Enables
Enable_SaveFigure = 1;

%% Plot
figure(1)

LineWidth = 1.5;
NumPoint = 500;
FigSize = [0.1 0.1 0.35 0.55];

set(gcf,'units','normalized','outerposition',FigSize);

% When lambda = 0, wb = (sqrt(2)-1)*60Hz = 24.8528Hz, which is the lower
% limit for wb.
F0 = 60;        % Fundamental frequency
fb(1) = (sqrt(2)-1)*F0;
fb(2) = 2/3*F0;
fb(3) = F0;
fb(4) = 1.5*F0;
fb(5) = 2*F0;

% FreqLower = min(fb);
% FreqUpper = max(fb);
FreqLower = 0.4*F0;
FreqUpper = 2*F0;

% ColorLower = [1,1,0]; ColorUpper = [1,0,0];       % yellow to red
% ColorLower = [1,0,1]; ColorUpper = [1,0,0];       % pink to red
ColorLower = [0,1,1]; ColorUpper = [0,0,1];       % light blue to dark blue
% ColorLower = [0,0,1]; ColorUpper = [0,0,0];       % blue to black

for k = 1:length(fb)
    if fb(k) <= F0
        AxisScale = 150;
    else
        AxisScale = 250;
    end
        
    x_vec = linspace(-AxisScale,0,NumPoint);
    y_vec = linspace(-AxisScale,AxisScale,NumPoint);
    y_vec = flip(y_vec);
    ConMat = GetContourMatrix(fb(k),x_vec,y_vec);
    
    FreqFactor = (fb(k)-FreqLower)/(FreqUpper-FreqLower);
    LineColor = (ColorUpper - ColorLower) * FreqFactor + ColorLower;
    
    contour(x_vec/F0,y_vec/F0,ConMat,[1,1],'LineWidth',LineWidth,'LineColor',LineColor); grid on; hold on;
end


xlabel('Real Part of $\lambda$ (pu)','interpreter','latex')
ylabel('Imaginary Part of $\lambda$ (pu)','interpreter','latex')

xlim([-3,0.5]);
xticks([-3,-2,-1,0,0.5])
ylim([-4,4]);
yticks([-4,-2,0,2,4]);

% Colorbar
ColorStep = 100;
GradRed     = linspace(ColorLower(1),ColorUpper(1),ColorStep)';
GradGreen   = linspace(ColorLower(2),ColorUpper(2),ColorStep)';
GradBlue    = linspace(ColorLower(3),ColorUpper(3),ColorStep)';
colormap([GradRed GradGreen GradBlue]);
colorbar();
caxis([FreqLower/F0,FreqUpper/F0]);

%% Plot pole
% Voltage node
if 1
% plot(0, 0, 'k.', 'markersize', 15); grid on; hold on;
scatter(0,0,'x','LineWidth',1,'MarkerEdgeColor',[0,0,0]); grid on; hold on;
end

% Current node
if 1
% pole_pu = load('H_pole_pu').pole_pu;
pole_pu = load('H_pole_pu_sweep').pole_pu;
for k = 1:length(pole_pu)
% plot(real(pole_pu(k)), imag(pole_pu(k)), 'k.', 'markersize', 15); grid on; hold on;
scatter(real(pole_pu{k}),imag(pole_pu{k}),'x','LineWidth',1,'MarkerEdgeColor',[0,0,0]); hold on; grid on;
end
end

%% Save
if Enable_SaveFigure
    print(gcf,'H_Contour.png','-dpng','-r600');
end
