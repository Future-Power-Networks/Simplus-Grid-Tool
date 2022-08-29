% Authos(s): Yitong Li

clear all
clc
close all

%%
load('SimIBR_Reference');
Reference = out.IbrData;
load('SimIBR_DC30Hz');
DC30Hz = out.IbrData;
load('SimIBR_AC150Hz');
AC150Hz = out.IbrData;
load('SimIBR_PLL30Hz');
PLL30Hz = out.IbrData;

time = Reference{1}.Values.Time;

idq_ref = Reference{1}.Values.Data;
w_ref_ = Reference{2}.Values.Data;
for i = 1:length(w_ref_)
    w_ref(i) = w_ref_(:,:,i);
end
vdc_ref_ = Reference{3}.Values.Data;
for i = 1:length(vdc_ref_)
    vdc_ref(i) = vdc_ref_(:,:,i);
end

idq_ac = AC150Hz{1}.Values.Data;
w_ac_ = AC150Hz{2}.Values.Data;
for i = 1:length(w_ac_)
    w_ac(i) = w_ac_(:,:,i);
end
vdc_ac_ = AC150Hz{3}.Values.Data;
for i = 1:length(vdc_ac_)
    vdc_ac(i) = vdc_ac_(:,:,i);
end

idq_dc = DC30Hz{1}.Values.Data;
w_dc_ = DC30Hz{2}.Values.Data;
for i = 1:length(w_dc_)
    w_dc(i) = w_dc_(:,:,i);
end
vdc_dc_ = DC30Hz{3}.Values.Data;
for i = 1:length(vdc_dc_)
    vdc_dc(i) = vdc_dc_(:,:,i);
end

idq_pll = PLL30Hz{1}.Values.Data;
w_pll_ = PLL30Hz{2}.Values.Data;
for i = 1:length(w_pll_)
    w_pll(i) = w_pll_(:,:,i);
end
vdc_pll_ = PLL30Hz{3}.Values.Data;
for i = 1:length(vdc_pll_)
    vdc_pll(i) = vdc_pll_(:,:,i);
end

%%
BLUE = [0, 0.4470, 0.7410];
RED = [0.8500, 0.3250, 0.0980];

FigN = 0;
LineWidth = 1;
EnableSave = 1;

time = time - 5.9;
XLim = [0,0.4];
YLim1 = [0.8,1.2];
YLim2 = [-1.5,1];
YLim3 = [0.93,1.07];
YTicks3 = [0.93,1,1.07];

FigureSize = [0.1 0.1 0.17 0.6];

FigN = FigN+1;
figure(FigN)
set(gcf,'units','normalized','outerposition',FigureSize);
subplot(3,1,1)
plot(time,w_ref,'-.','LineWidth',LineWidth,'color',[.5 .5 .5]); grid on; hold on;
ylabel('$\omega_8$ (pu)','interpreter','latex')
xlim(XLim)
ylim(YLim1)
subplot(3,1,2)
plot(time,idq_ref(:,1),'-.','LineWidth',LineWidth,'color',[.5 .5 .5]); grid on; hold on;
plot(time,idq_ref(:,2),'-.','LineWidth',LineWidth,'color',[.5 .5 .5]); grid on; hold on;
ylabel('$i_{dq8}$ (pu)','interpreter','latex')
xlim(XLim)
ylim(YLim2)
subplot(3,1,3)
plot(time,vdc_ref,'-.','LineWidth',LineWidth,'color',[.5 .5 .5]); grid on; hold on;
ylabel('$v_{dc8}$ (pu)','interpreter','latex')
xlabel('Times (s)','interpreter','latex')
xlim(XLim)
ylim(YLim3)
yticks(YTicks3)
if EnableSave
    print(gcf,'Simulation/Figure/SimIBR_Ref.png','-dpng','-r600');
end

FigN = FigN+1;
figure(FigN)
set(gcf,'units','normalized','outerposition',FigureSize);
subplot(3,1,1)
plot(time,w_ref,'-.','LineWidth',LineWidth,'color', [.5 .5 .5]); grid on; hold on;
plot(time,w_ac,'LineWidth',LineWidth,'color',BLUE); grid on; hold on;
ylabel('$\omega_8$ (pu)','interpreter','latex')
xlim(XLim)
ylim(YLim1)
subplot(3,1,2)
plot(time,idq_ref(:,1),'-.','LineWidth',LineWidth,'color',[.5 .5 .5]); grid on; hold on;
plot(time,idq_ref(:,2),'-.','LineWidth',LineWidth,'color',[.5 .5 .5]); grid on; hold on;
plot(time,idq_ac(:,1),'LineWidth',LineWidth,'color',BLUE); grid on; hold on;
plot(time,idq_ac(:,2),'LineWidth',LineWidth,'color',RED); grid on; hold on;
ylabel('$i_{dq8}$ (pu)','interpreter','latex')
xlim(XLim)
ylim(YLim2)
subplot(3,1,3)
plot(time,vdc_ref,'-.','LineWidth',LineWidth,'color',[.5 .5 .5]); grid on; hold on;
plot(time,vdc_ac,'-.','LineWidth',LineWidth,'color',BLUE); grid on; hold on;
ylabel('$v_{dc8}$ (pu)','interpreter','latex')
xlabel('Times (s)','interpreter','latex')
xlim(XLim)
ylim(YLim3)
yticks(YTicks3)
if EnableSave
    print(gcf,'Simulation/Figure/SimIBR_AC150Hz.png','-dpng','-r600');
end

FigN = FigN+1;
figure(FigN)
set(gcf,'units','normalized','outerposition',FigureSize);
subplot(3,1,1)
plot(time,w_ref,'-.','LineWidth',LineWidth,'color', [.5 .5 .5]); grid on; hold on;
plot(time,w_pll,'LineWidth',LineWidth,'color',BLUE); grid on; hold on;
ylabel('$\omega_8$ (pu)','interpreter','latex')
xlim(XLim)
ylim(YLim1)
subplot(3,1,2)
plot(time,idq_ref(:,1),'-.','LineWidth',LineWidth,'color',[.5 .5 .5]); grid on; hold on;
plot(time,idq_ref(:,2),'-.','LineWidth',LineWidth,'color',[.5 .5 .5]); grid on; hold on;
plot(time,idq_pll(:,1),'LineWidth',LineWidth,'color',BLUE); grid on; hold on;
plot(time,idq_pll(:,2),'LineWidth',LineWidth,'color',RED); grid on; hold on;
ylabel('$i_{dq8}$ (pu)','interpreter','latex')
xlim(XLim)
ylim(YLim2)
subplot(3,1,3)
plot(time,vdc_ref,'-.','LineWidth',LineWidth,'color',[.5 .5 .5]); grid on; hold on;
plot(time,vdc_pll,'-.','LineWidth',LineWidth,'color',BLUE); grid on; hold on;
ylabel('$v_{dc8}$ (pu)','interpreter','latex')
xlabel('Times (s)','interpreter','latex')
xlim(XLim)
ylim(YLim3)
yticks(YTicks3)
if EnableSave
    print(gcf,'Simulation/Figure/SimIBR_PLL30Hz.png','-dpng','-r600');
end

FigN = FigN+1;
figure(FigN)
set(gcf,'units','normalized','outerposition',FigureSize);
subplot(3,1,1)
plot(time,w_ref,'-.','LineWidth',LineWidth,'color', [.5 .5 .5]); grid on; hold on;
plot(time,w_dc,'LineWidth',LineWidth,'color',BLUE); grid on; hold on;
ylabel('$\omega_8$ (pu)','interpreter','latex')
xlim(XLim)
ylim(YLim1)
subplot(3,1,2)
plot(time,idq_ref(:,1),'-.','LineWidth',LineWidth,'color',[.5 .5 .5]); grid on; hold on;
plot(time,idq_ref(:,2),'-.','LineWidth',LineWidth,'color',[.5 .5 .5]); grid on; hold on;
plot(time,idq_dc(:,1),'LineWidth',LineWidth,'color',BLUE); grid on; hold on;
plot(time,idq_dc(:,2),'LineWidth',LineWidth,'color',RED); grid on; hold on;
ylabel('$i_{dq8}$ (pu)','interpreter','latex')
xlim(XLim)
ylim(YLim2)
subplot(3,1,3)
plot(time,vdc_ref,'-.','LineWidth',LineWidth,'color',[.5 .5 .5]); grid on; hold on;
plot(time,vdc_dc,'-.','LineWidth',LineWidth,'color',BLUE); grid on; hold on;
ylabel('$v_{dc8}$ (pu)','interpreter','latex')
xlabel('Times (s)','interpreter','latex')
xlim(XLim)
ylim(YLim3)
yticks(YTicks3)
if EnableSave
    print(gcf,'Simulation/Figure/SimIBR_DC30Hz.png','-dpng','-r600');
end

