% Author(s): Yitong Li

clc
close all

FigN = 0;

FigN = FigN+1;
figure(FigN);
[PfAllValue,PfAllState,Mode] = ModalAnalysisDss(ObjGsysSs,594,30,FigN,0);
clear('PfValue')
clear('PfState')
pie(PfAllValue,PfAllState);

FigN = FigN+1;
figure(FigN);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.25 0.4]);
PfValue(1) = sum(PfAllValue([1:12]));
PfValue(2) = sum(PfAllValue([13:end]));
PfState{1} = 'IBR10';
PfState{2} = 'Others';
pie(PfValue,PfState);

if 1
    print(gcf,'Simulation118Bus/118Bus_Participation.png','-dpng','-r600');
end