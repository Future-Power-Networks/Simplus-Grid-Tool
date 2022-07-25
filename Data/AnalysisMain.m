% Author(s): Yitong Li

clear all
clc
close all

GflYsym = load('Data/GflYsym');
GfmYsym = load('Data/GfmYsym');

GflYss = load('Data/GflYss');
GfmYss = load('Data/GfmYss');

% GflYsym = load('Data/GflDcLinkYsym');
% GfmYsym = load('Data/SgYsym');

OmegaPositive = logspace(-1,4,1e3)*2*pi;


%%
figure(1)

set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.8]);

SimplusGT.bode_c(GflYsym.Ysym{2}(1,1),1j*OmegaPositive,'PhaseOn',1); hold on;
SimplusGT.bode_c(GfmYsym.Ysym{2}(1,1),1j*OmegaPositive,'PhaseOn',1,'PhaseShift',2*pi); hold on;

xlabel('Frequency (Hz)','interpreter','Latex')

subplot(2,1,1)
yticks([1e-3,1e-2,1e-1,1e0,1e1,1e2])
ylim([1e-3,1e2])
ylabel('Admittance (pu)','interpreter','Latex')
legend('Grid-Following','Grid-Forming')
subplot(2,1,2)
ylabel('Phase (Degree)','interpreter','Latex')

if 0
figure(1)
print(gcf,'Data/AdmittanceSpectrum.png','-dpng','-r600');
end

%%
figure(2)

set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.8]);

GflZss = inv(GflYss.Yss{2});
GfmZss = inv(GfmYss.Yss{2});

options = bodeoptions;
options.FreqUnits = 'Hz';

a1 = bode(GflZss(1,1),options); grid on; hold on;
a2 = bode(GfmZss(1,1),options); grid on; hold on;

if 0
figure(2)
print(gcf,'Data/ImpedanceSpectrum.png','-dpng','-r600');
end


% SimplusGT.bode_c(GfmZsym(1,1),1j*OmegaPositive,'PhaseOn',1); hold on;