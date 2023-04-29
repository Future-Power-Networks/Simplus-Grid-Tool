% Author(s): Yitong Li

clear all
clc
close all

EnableSave = 1;

%%
% Load data
EigVecGFM{1} = load('EigVec_GFM_P0_Q0').EigVec;
EigVecGFM{2} = load('EigVec_GFM_P05_Q0').EigVec;
EigVecGFM{3} = load('EigVec_GFM_Pn05_Q0').EigVec;
EigVecGFM{4} = load('EigVec_GFM_P0_Q05').EigVec;
EigVecGFM{5} = load('EigVec_GFM_P0_Qn05').EigVec;

EigVecGFL{1} = load('EigVec_GFL_P0_Q0').EigVec;
EigVecGFL{2} = load('EigVec_GFL_P05_Q0').EigVec;
EigVecGFL{3} = load('EigVec_GFL_Pn05_Q0').EigVec;
EigVecGFL{4} = load('EigVec_GFL_P0_Q05').EigVec;
EigVecGFL{5} = load('EigVec_GFL_P0_Qn05').EigVec;

%%
% Convert rad/s to Hz
for i = 1:length(EigVecGFM)
    EigVecGFM{i} = EigVecGFM{i}/2/pi;
end

for i = 1:length(EigVecGFL)
    EigVecGFL{i} = EigVecGFL{i}/2/pi;
end

%%
% Plot
figure(1)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.5]);
for i = 1:length(EigVecGFM)
    PlotPole(EigVecGFM{i});
end
legend('P=0,Q=0','P=0.5,Q=0','P=-0.5,Q=0','P=0,Q=0.5','P=0,Q=-0.5','location','west')
axis([-2,2,-20,20])
if EnableSave
    print(gcf,'SimOperatingPoint/PoleGfm.png','-dpng','-r600');
end

figure(2)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.5]);
for i = 1:length(EigVecGFL)
    PlotPole(EigVecGFL{i});
end
legend('P=0,Q=0','P=0.5,Q=0','P=-0.5,Q=0','P=0,Q=0.5','P=0,Q=-0.5','location','west')
axis([-4,4,-30,30])
if EnableSave
    print(gcf,'SimOperatingPoint/PoleGfl.png','-dpng','-r600');
end

%%
function PlotPole(Eigvalue)
    scatter(real(Eigvalue),imag(Eigvalue),'x','LineWidth',1.5); hold on; grid on;
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
end