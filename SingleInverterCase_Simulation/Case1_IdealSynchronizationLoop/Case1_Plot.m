% Author(s): Yitong Li

%%
clear all
clc
close all

%%
% Load data
Case1_FDC = load('Case1_FDC').Case1_FDC;
Case1_PLL = load('Case1_PLL').Case1_PLL;

% time
time = Case1_FDC.time;

% Organize data
ig_FDC = Case1_FDC.signals(1).values;
idq_FDC = Case1_FDC.signals(2).values;
w_FDC = Case1_FDC.signals(3).values;
itheta_FDC = Case1_FDC.signals(5).values;

[~,c_FDC] = size(idq_FDC);

for i = 1:c_FDC/2
    id_FDC(:,i) = idq_FDC(:,(i*2-1));
    iq_FDC(:,i) = idq_FDC(:,(i*2));
end

vg_PLL = Case1_PLL.signals(1).values;
vdq_PLL = Case1_PLL.signals(2).values;
w_PLL = Case1_PLL.signals(3).values;
vtheta_PLL = Case1_PLL.signals(5).values;

[~,c_PLL] = size(vdq_PLL);

for i = 1:c_PLL/2
    vd_PLL(:,i) = vdq_PLL(:,(i*2-1));
    vq_PLL(:,i) = vdq_PLL(:,(i*2));
end

%%
% Plot

fn = 0;
LineWidth = 1;
ylim_vi = [-1.2,1.2];
yticks_vi = [-1,0,1];
enable_save = 1;

fn = fn+1;
figure(fn)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.2 0.5]); % position of left-bottow cornor + length/depth of figure
subplot(4,1,1)
plot(time,id_FDC,'LineWidth',LineWidth); grid on;
ylabel('$i_d$ (pu)','interpreter','latex')
ylim(ylim_vi);
yticks(yticks_vi);
subplot(4,1,2)
plot(time,iq_FDC,'LineWidth',LineWidth); grid on;
ylabel('$i_q$ (pu)','interpreter','latex')
ylim(ylim_vi);
yticks(yticks_vi);
subplot(4,1,3)
plot(time,w_FDC,'LineWidth',LineWidth); grid on;
ylabel('$\omega$ (pu)','interpreter','latex')
subplot(4,1,4)
plot(time,itheta_FDC,'LineWidth',LineWidth); grid on;
ylabel('$\theta_i$ ($^o$)','interpreter','latex')
ylim([-270,90]);
yticks([-270,-180,-90,0,90]);

% legend('$\theta_{i0}=45^o$','$\theta_{i0}=0^o$','$\theta_{i0}=-45^o$','$\theta_{i0}=-90^o$','$\theta_{i0}=-135^o$','$\theta_{i0}=-180^o$','$\theta_{i0}=-225^o$','interpreter','latex')
xlabel('Time (s)')
if enable_save
    print(gcf,'Case1_SimFDC.png','-dpng','-r600'); 
end

fn = fn+1;
figure(fn)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.2 0.5]); % position of left-bottow cornor + length/depth of figure
subplot(4,1,1)
plot(time,vd_PLL,'LineWidth',LineWidth); grid on;
ylabel('$v_d$ (pu)','interpreter','latex')
ylim(ylim_vi);
yticks(yticks_vi);
subplot(4,1,2)
plot(time,vq_PLL,'LineWidth',LineWidth); grid on;
ylabel('$v_q$ (pu)','interpreter','latex')
ylim(ylim_vi);
yticks(yticks_vi);
subplot(4,1,3)
plot(time,w_PLL,'LineWidth',LineWidth); grid on;
ylabel('$\omega$ (pu)','interpreter','latex')
subplot(4,1,4)
plot(time,vtheta_PLL,'LineWidth',LineWidth); grid on;
ylabel('$\theta_v$ ($^o$)','interpreter','latex')
ylim([-180,180]);
yticks([-180,-90,0,90,180]);
xlabel('Time (s)','interpreter','latex')

if enable_save
    print(gcf,'Case1_SimPLL.png','-dpng','-r600'); 
end