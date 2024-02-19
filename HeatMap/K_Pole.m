% This function analyzes the pole of 68 bus system

clear all
clc
close all

mfile_name = mfilename('fullpath');
[RootPath,~,~]  = fileparts(mfile_name);
cd(RootPath);

%% Enables
Enable_SaveFigure = 1;

%% Load data
pole{1}    = load('K_SG_Load_pole_T12cl').pole_T12cl;
pole{2}    = load('K_SG_IBR_Load_pole_T12cl').pole_T12cl;
pole{3}    = load('K_SG_IBR_pole_T12cl').pole_T12cl;
pole{4}    = load('K_SG_IBR_17_pole_T12cl').pole_T12cl;

%% Plot settings
FigNum = 0;
LineWidth = 1.5;

y_Limit = [-4,4];
y_Ticks = [-4,-2,0,2,4];
x_Limit = [-8,2];
x_Ticks = [-8,-6,-4,-2,0,2];

FigSize = [0.1 0.1 0.35 0.5];

%% Plot
for i = 1:4
    FigNum = FigNum + 1;
    figure(FigNum)
    set(gcf,'units','normalized','outerposition',FigSize);
    scatter(real(pole{i}),imag(pole{i}),'x','LineWidth',LineWidth); hold on; grid on;
    xlabel('Real Part (Hz)','interpreter','latex')
    ylabel('Imagary Part (Hz)','interpreter','latex')
    if i ~= 1
    xlim(x_Limit);
    xticks(x_Ticks)
    end
    ylim(y_Limit);
    yticks(y_Ticks);
    if Enable_SaveFigure == 1
        print(gcf,['K_Pole_' num2str(i) '.png'],'-dpng','-r600');
    end
end