% This function does the theoratical analysis of
% single-machine-infinite-bus case.

% Author(s): Yitong Li

% Notes:
% 
% Matlab defaultly has 7 different colors. The 1st one is dark blue,
% and the 7th one is dark red.

%%
clear all
clc
close all

%%
s = sym('s');       % Laplace operator

Wbase = 2*pi*50;    % (rad/s)
Vbase = 1;
Sbase = 1;
Ibase = Sbase/Vbase;
Zbase = Vbase/Ibase;
Ybase = 1/Zbase;

omega_p = logspace(0,3,500)*2*pi;
omega_pn = [-flip(omega_p),omega_p];

enable_save = 0;
lw = 1;

%% Droop grid-forming inverter
fn = 100;

% ### Plot bandwidth
if 0
[~,~,~,~,~,~,~,Gv_cl,Gi_cl] = CalcSwingDroop();
Gv_cl = SimplexPS.tf2sym(Gv_cl);
Gi_cl = SimplexPS.tf2sym(Gi_cl);
fn = fn+1; 
figure(fn)
SimplexPS.bode_c(Gv_cl,1j*omega_pn,2*pi,'PhaseOn',0); hold on;
SimplexPS.bode_c(Gi_cl,1j*omega_pn,2*pi,'PhaseOn',0); hold on;
end

% ### Change series impedance R+jX
if 0
fn = fn+1; 
figure(fn)
clear('root','root_prime','XRratio','Zline','Xline','Rline');

% XRratio = 5;
% Zline = linspace(0.2,0.1,7);
% [Xline,Rline] = CalcLine(Zline,XRratio);
Xline = linspace(0.3,0.1,7);
Rline = Xline/5;

set(gcf,'units','normalized','outerposition',[0.1 0.1 0.18 0.4]);
for i = 1:length(Xline)
    [~,~,~,~,root{i},root_prime{i}] = CalcSwingDroop('Xline',Xline(i),'Rline',Rline(i));
    
    subplot(2,1,1)
    scatter(real(root_prime{i}),imag(root_prime{i}),'x','LineWidth',lw); hold on; grid on;
  	xlabel('Real (Hz)','interpreter','latex')
    ylabel('Imaginary (Hz)','interpreter','latex')
    xlim([-200,50]);
    xticks([-200,-150,-100,-50,0,50])
    ylim([-2500,2500]);
    yticks([-2500,-1250,0,1250,2500]);
    % plot([0,0],[-2500,2500],'--k'); hold on; grid on;
    
    subplot(2,1,2)
    scatter(real(root_prime{i}),imag(root_prime{i}),'x','LineWidth',lw); hold on; grid on;
 	xlabel('Real (Hz)','interpreter','latex')
    ylabel('Imaginary (Hz)','interpreter','latex')
    title('Zoomed-in Plot','interpreter','latex')
	xlim([-20,10]);
    ylim([-40,40]);
    yticks([-40,-20,0,20,40]);
end
if enable_save
    print(gcf,'Case2_GFM_Z.png','-dpng','-r600');
end

end

% ### Change droop gain
if 0
fn = fn+1; 
figure(fn)
clear('root','root_prime','XRratio','Zline','Xline','Rline');
 
m = 0.05*linspace(1,4,7);

set(gcf,'units','normalized','outerposition',[0.1 0.1 0.18 0.4]);
for i = 1:length(m)
    [~,~,~,~,root{i},root_prime{i}] = CalcSwingDroop('m',m(i));
    
    subplot(2,1,1)
    scatter(real(root_prime{i}),imag(root_prime{i}),'x','LineWidth',lw); hold on; grid on;
    xlabel('Real (Hz)','interpreter','latex')
    ylabel('Imaginary (Hz)','interpreter','latex')
    xlim([-200,50]);
    xticks([-200,-150,-100,-50,0,50])
    ylim([-2500,2500]);
    yticks([-2500,-1250,0,1250,2500]);
    
    subplot(2,1,2)
    scatter(real(root_prime{i}),imag(root_prime{i}),'x','LineWidth',lw); hold on; grid on;
    xlabel('Real (Hz)','interpreter','latex')
    ylabel('Imaginary (Hz)','interpreter','latex')
    title('Zoomed-in Plot','interpreter','latex')
	xlim([-20,10]);
    ylim([-40,40]);
    yticks([-40,-20,0,20,40]);
end
if enable_save
    print(gcf,'Case2_GFM_m.png','-dpng','-r600');
end
end

% ### Change voltage loop gain
if 0
fn = fn+1; 
figure(fn)
clear('root','root_prime','XRratio','Zline','Xline','Rline');

w_v = 250*2*pi*linspace(1,3/5,7);

set(gcf,'units','normalized','outerposition',[0.1 0.1 0.18 0.4]);
for i = 1:length(w_v)
    [~,~,~,~,root{i},root_prime{i}] = CalcSwingDroop('w_v',w_v(i));
    
    subplot(2,1,1)
    scatter(real(root_prime{i}),imag(root_prime{i}),'x','LineWidth',lw); hold on; grid on;
    xlabel('Real (Hz)','interpreter','latex')
    ylabel('Imaginary (Hz)','interpreter','latex')
    xlim([-200,50]);
    xticks([-200,-150,-100,-50,0,50])
    ylim([-2500,2500]);
    yticks([-2500,-1250,0,1250,2500]);
    
    subplot(2,1,2)
    scatter(real(root_prime{i}),imag(root_prime{i}),'x','LineWidth',lw); hold on; grid on;
    xlabel('Real (Hz)','interpreter','latex')
    ylabel('Imaginary (Hz)','interpreter','latex')
    title('Zoomed-in Plot','interpreter','latex')
	xlim([-20,10]);
    ylim([-40,40]);
    yticks([-40,-20,0,20,40]);
end
if enable_save
    print(gcf,'Case2_GFM_w_v.png','-dpng','-r600');
end
end

%% PLL grid-following inverter
% Initialize figure index
fn = 200;

% ### Plot bandwidth
if 0
[~,~,~,~,~,~,~,Gi_cl] = CalcSwingPLL();
Gi_cl = simplify(Gi_cl);
Gi_cl = SimplexPS.tf2sym(Gi_cl);
fn = fn+1; 
figure(fn)
SimplexPS.bode_c(Gi_cl,1j*omega_pn,2*pi,'PhaseOn',0);
end

% ### Change line impedance
if 1
    
fn = fn+1; 
figure(fn)
clear('root','root_prime','XRratio','Zline','Xline','Rline');     

Xline = linspace(0.4,0.8,10);
Rline = Xline/10;

% set(gcf,'units','normalized','outerposition',[0.1 0.1 0.18 0.4]);
for i = 1:length(Xline)
	[~,~,~,~,~,root_prime{i},pole_sys{i}] = CalcSwingPLL('Xline',Xline(i),'Rline',Rline(i));
    
    subplot(2,1,1)
    scatter(real(pole_sys{i}),imag(pole_sys{i}),'x','LineWidth',lw); hold on; grid on;
    xlabel('Real (Hz)','interpreter','latex')
    ylabel('Imaginary (Hz)','interpreter','latex')
 	xlim([-200,50]);
    xticks([-200,-150,-100,-50,0,50])
    ylim([-2500,2500]);
    yticks([-2500,-1250,0,1250,2500]);
    
    subplot(2,1,2)
   	scatter(real(pole_sys{i}),imag(pole_sys{i}),'x','LineWidth',lw); hold on; grid on;
    % scatter(real(pole_sys{i}),imag(pole_sys{i}),'x','LineWidth',lw); hold on; grid on;
    xlabel('Real (Hz)','interpreter','latex')
    ylabel('Imaginary (Hz)','interpreter','latex')
    title('Zoomed-in Plot','interpreter','latex')
   	xlim([-10,5]);
    ylim([-60,60]);
    % yticks([-40,-20,0,20,40]);
end


if enable_save
    print(gcf,'Case2_GFF_Y.png','-dpng','-r600');
end
end

% ### Change PLL bandwidth
if 0
    
fn = fn+1; 
figure(fn)
clear('root','root_prime','XRratio','Zline','Xline','Rline');    

w_pll = 15*2*pi * linspace(1,4,7);

set(gcf,'units','normalized','outerposition',[0.1 0.1 0.18 0.4]);
for i = 1:length(w_pll)
	[~,~,~,~,~,root_prime{i,1},pole_sys{i,1}] = CalcSwingPLL('w_pll',w_pll(i));
    
    subplot(2,1,1)
    scatter(real(root_prime{i,1}),imag(root_prime{i,1}),'x','LineWidth',lw); hold on; grid on;
    xlabel('Real (Hz)','interpreter','latex')
    ylabel('Imaginary (Hz)','interpreter','latex')
  	xlim([-200,50]);
    xticks([-200,-150,-100,-50,0,50])
    ylim([-2500,2500]);
    yticks([-2500,-1250,0,1250,2500]);
    
    subplot(2,1,2)
    scatter(real(root_prime{i,1}),imag(root_prime{i,1}),'x','LineWidth',lw); hold on; grid on;
    xlabel('Real (Hz)','interpreter','latex')
    ylabel('Imaginary (Hz)','interpreter','latex')
    title('Zoomed-in Plot','interpreter','latex')
  	xlim([-10,5]);
    ylim([-40,40]);
    yticks([-40,-20,0,20,40]);
end


if enable_save
    print(gcf,'Case2_GFF_w_pll.png','-dpng','-r600');
end
end

% ### Change current loop bandwidth
if 0
    
fn = fn+1; 
figure(fn)
clear('root','root_prime','XRratio','Zline','Xline','Rline');    

w_i = 250*2*pi * linspace(1,3/5,7);

set(gcf,'units','normalized','outerposition',[0.1 0.1 0.18 0.4]);
for i = 1:length(w_i)
	[~,~,~,~,~,root_prime{i,1},pole_sys{i,1}] = CalcSwingPLL('w_i',w_i(i));
    
    subplot(2,1,1)
    scatter(real(root_prime{i,1}),imag(root_prime{i,1}),'x','LineWidth',lw); hold on; grid on;
    xlabel('Real (Hz)','interpreter','latex')
    ylabel('Imaginary (Hz)','interpreter','latex')
  	xlim([-200,50]);
    xticks([-200,-150,-100,-50,0,50])
    ylim([-2500,2500]);
    yticks([-2500,-1250,0,1250,2500]);
    
    subplot(2,1,2)
    scatter(real(root_prime{i,1}),imag(root_prime{i,1}),'x','LineWidth',lw); hold on; grid on;
    xlabel('Real (Hz)','interpreter','latex')
    ylabel('Imaginary (Hz)','interpreter','latex')
    title('Zoomed-in Plot','interpreter','latex')
  	xlim([-10,5]);
    ylim([-40,40]);
    yticks([-40,-20,0,20,40]);
end


if enable_save
    print(gcf,'Case2_GFF_w_i.png','-dpng','-r600');
end
end