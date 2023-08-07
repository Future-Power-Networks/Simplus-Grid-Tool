% This script plots the gamma and damping analysis

clear all
clc
close all

mfile_name = mfilename('fullpath');
[RootPath,~,~]  = fileparts(mfile_name);
cd(RootPath);

ColorRGB();

%% Enables
Enable_SaveFigure = 1;

%% Load data
GammaSimData_Unstable = load('Gamma_SimData_SG_Unstable').Gamma_SimData_SG;
GammaSimData_Stable = load('Gamma_SimData_SG_Stable').Gamma_SimData_SG;

% Unstable case
Time_us = GammaSimData_Unstable.time;
w_us__ = GammaSimData_Unstable.signals(5).values;
for i = 1:length(w_us__)
    w_us(i) = w_us__(:,:,i);
end
theta_us__ = GammaSimData_Unstable.signals(6).values;
for i = 1:length(theta_us__)
    theta_us(i) = theta_us__(:,:,i);
end

% Stable case
Time_st = GammaSimData_Stable.time;
w_st__ = GammaSimData_Stable.signals(5).values;
for i = 1:length(w_st__)
    w_st(i) = w_st__(:,:,i);
end
theta_st__ = GammaSimData_Stable.signals(6).values;
for i = 1:length(theta_st__)
    theta_st(i) = theta_st__(:,:,i);
end

%% Calc
% Unstable case
TimeRef_us = find(Time_us == 3);
TimeStart_us = find(Time_us == 3.35);
TimeEnd_us = find(Time_us == 3.5);
w_ref_us = w_us(TimeRef_us-100);
theta_ref_us = theta_us(TimeRef_us-100);
w_us = w_us(TimeStart_us:TimeEnd_us) - w_ref_us;                     % The first few points with "noises" are discarded.
theta_us = theta_us(TimeStart_us:TimeEnd_us) - theta_ref_us;

% Normalize
if 0                        % We do not normalize the results by default
    W0 = 2*pi*60;
    J = 0.05*2/W0;
    K = 2.8253;            % Measured by runing MainSyn
    if 0
        w_us = w_us*W0*sqrt(J);
        theta_us = theta_us*sqrt(K);
    else
        w_us = w_us;
        theta_us = theta_us*sqrt(K)/sqrt(J)/W0;
    end
end

% Stable case
TimeRef_st = find(Time_st == 3);
TimeStart_st = find(Time_st == 3.35);
TimeEnd_st = find(Time_st == 3.5);
w_ref_st = w_st(TimeRef_st-100);
theta_ref_st = theta_st(TimeRef_st-100);
w_st = w_st(TimeStart_st:TimeEnd_st) - w_ref_st;                     % The first few points with "noises" are discarded.
theta_st = theta_st(TimeStart_st:TimeEnd_st) - theta_ref_st;

%% Plot settings
FigNum = 0;
LineWidth = 1.2;

FigRowMax = 1;
FigColumnMax = 2;

x_Limit = [-0.04,0.04];
x_Ticks = [-0.04,-0.02,0,0.02,0.04];
y_Limit = [-0.12,0.12];
y_Ticks = [-0.12,-0.06,0,0.06,0.12];

%% Plot

FigNum = FigNum + 1;
figure(FigNum)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.35 0.63]);

% Unstable phase portrait
h1 = plot(w_us,theta_us,'LineWidth',LineWidth,'Color',RgbBlue); grid on; hold on;

% Stable case
h2 = plot(w_st,theta_st,'LineWidth',LineWidth,'Color',RgbRed); grid on; hold on;

% Legend
% legend({'Unstable','Stable'},'interpreter','latex')
legend({'Unstable','Stable'})

% Settings
xlabel('$\Delta \omega$ (pu)','interpreter','latex')
ylabel('$\Delta \theta$ (rad)','interpreter','latex')
xlim(x_Limit);
xticks(x_Ticks)
ylim(y_Limit);
yticks(y_Ticks);

% Unstable arrows
StepNum = 100;
StepSize = floor(length(w_us)/StepNum);
for i = 1:(StepNum-1)
    w1 = w_us(i*StepSize+1);
    theta1 = theta_us(i*StepSize+1);
    ScaleFactor(i) = abs(w1)^2*370;
    h_us_quiver = quiver(w1,theta1,w1,theta1,ScaleFactor(i),'Color',RgbBlue,'MaxHeadSize',0.8); grid on; hold on;
    h_us_quiver.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

% Stable arrows
StepNum = 80;
StepSize = floor(length(w_st)/StepNum);
for i = 1:(StepNum-1)
    w1 = w_st(i*StepSize+1);
    theta1 = theta_st(i*StepSize+1);
    ScaleFactor(i) = abs(w1)^2*1800;
    h_st_quiver = quiver(w1,theta1,-w1,-theta1,ScaleFactor(i),'Color',RgbRed,'MaxHeadSize',0.8); grid on; hold on;
    h_st_quiver.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

if Enable_SaveFigure
    print(gcf,'Gamma_Sim_SG.png','-dpng','-r600');
end



