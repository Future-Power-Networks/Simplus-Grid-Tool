% Author(s): Yitong Li

clc
close all

FigN = 0;

FigN = FigN+1;
figure(FigN);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.4]);
[PfAllValue,PfAllState,Mode] = ModalAnalysisDss(ObjGsysSs,594,30,FigN,0);
clear('PfValue')
clear('PfState')
pie(PfAllValue,PfAllState);

% PfValue(1) = sum(PfAllValue([14:22,27]));
% PfState{1} = 'IBR2';
% PfValue(2) = sum(PfAllValue([1:6,12,24:26,29,30]));
% PfState{2} = 'IBR7';
% PfValue(3) = sum(PfAllValue([7:11,13,23,28]));
% PfState{3} = 'IBR10';
% PfValue(4) = sum(PfAllValue(31:end));
% PfState{4} = 'Others';
% 
% subplot(1,2,1)
% % title('Mode -0.44+j36.79 Hz')
% % pie(PfAllValue,PfAllState);
% pie(PfValue,PfState);
% 
% [PfAllValue,PfAllState,Mode] = ModalAnalysisDss(ObjGsysSs,243,25,FigN,0);
% clear('PfValue')
% clear('PfState')
% PfValue(1) = sum(PfAllValue([7:16,19,20,22:24]));
% PfState{1} = 'IBR2';
% PfValue(2) = sum(PfAllValue([1:6,17,18,21,25]));
% PfState{2} = 'IBR10';
% PfValue(3) = sum(PfAllValue(26:end));
% PfState{3} = 'Others';
% 
% subplot(1,2,2)
% % title('Mode -0.15+j39.13 Hz')
% % pie(PfAllValue,PfAllState);
% pie(PfValue,PfState);

if 0
    print(gcf,'Simulation118Bus/118Bus_Participation.png','-dpng','-r600');
end