% This function analyzes the participation

clear all
clc
close all

mfile_name = mfilename('fullpath');
[RootPath,~,~]  = fileparts(mfile_name);
cd(RootPath);

%% Enables
Enable_SaveFigure = 1;

%% Load data
KH    = load('K_SG_IBR_KH').KH;
Order = load('K_SG_IBR_Order').Order_New_VInoF;

%% Calc
[phi,xi,psi] = eig(KH);
phi_inv = inv(phi);
xi_diag = diag(xi);
[~,xi_min_index] = min(real(xi_diag));
xi_min = xi_diag(xi_min_index);

phi_r_xi = phi(:,xi_min_index);
phi_l_xi = transpose(phi_inv(xi_min_index,:));
FiedlerVec = phi_r_xi.*phi_l_xi;
FiedlerAbs = abs(FiedlerVec);

%% Figure Settings
Index17 = find(Order == 17);
FigSize = [0.1 0.1 0.3 0.4];

%% Plot
figure(1)
set(gcf,'units','normalized','outerposition',FigSize);
h = bar(Order,FiedlerAbs,'FaceColor',[0, 0.4470, 0.7410]); grid on; hold on;
bar(Order(Index17),FiedlerAbs(Index17),'FaceColor',[0.8500, 0.3250, 0.0980]); grid on; hold on;
ylim([0,0.8]);
yticks([0,0.2,0.4,0.6,0.8]);
xlabel('Bus Number','interpreter','latex')
ylabel('Participation Factor','interpreter','latex')

if Enable_SaveFigure == 1
    print(gcf,['K_Bar.png'],'-dpng','-r600');
end