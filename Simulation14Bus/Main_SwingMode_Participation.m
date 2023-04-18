% Author(s): Yitong Li

clc
close all

EnableSave = 1;

FigN = 0;

FigN = FigN+1;
figure(FigN);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.2 0.4]);
[PfAllValue,PfStateStr,Mode] = ModalAnalysisDss(ObjGsysSs,72,6,FigN,0);
for i = 1:length(PfStateStr)
    PfStateStr{i} = strrep(PfStateStr{i},'1','_1');
    PfStateStr{i} = strrep(PfStateStr{i},'3','_3');
    PfStateStr{i} = strrep(PfStateStr{i},'6','_6');
end
PfStateStr{3} = '';
PfStateStr{4} = '';
PfStateStr{5} = '';
PfStateStr{6} = '';
pie(PfAllValue,PfStateStr);
if EnableSave
    print(gcf,'Simulation/Figure/SwingModeHF.png','-dpng','-r600');
end

FigN = FigN+1;
figure(FigN);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.2 0.4]);
[PfAllValue,PfStateStr,Mode] = ModalAnalysisDss(ObjGsysSs,74,6,FigN,0);
for i = 1:length(PfStateStr)
    PfStateStr{i} = strrep(PfStateStr{i},'1','_1');
    PfStateStr{i} = strrep(PfStateStr{i},'3','_3');
    PfStateStr{i} = strrep(PfStateStr{i},'6','_6');
end
PfStateStr{5} = '';
PfStateStr{6} = '';
pie(PfAllValue,PfStateStr);
if EnableSave
    print(gcf,'Simulation/Figure/SwingModeLF.png','-dpng','-r600');
end