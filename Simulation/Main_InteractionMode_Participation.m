% Author(s): Yitong Li

clc
close all

EnableSave = 1;

FigN = 0;

FigN = FigN+1;
figure(FigN);
[PfAllValue,PfStateStr,~] = ModalAnalysisDss(ObjGsysSs,82,8,FigN,0);
pie(PfAllValue,PfStateStr);

FigN = FigN+1;
figure(FigN);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.2 0.4]);
[PfAllValue,PfStateStr,~] = ModalAnalysisDss(ObjGsysSs,82,8,FigN,0);

% Update value and name
PfAllValueNew(1) = PfAllValue(1);
PfStateStrNew{1} = 'i_{d8,ki}';

PfAllValueNew(2) = PfAllValue(2);
PfStateStrNew{2} = 'i_{q8,ki}';

PfAllValueNew(3) = PfAllValue(3)+PfAllValue(4);
PfStateStrNew{3} = 'i_{d8}';

PfAllValueNew(4) = PfAllValue(5)+PfAllValue(6);
PfStateStrNew{4} = 'i_{q8}';

PfAllValueNew(5) = PfAllValue(7);
PfStateStrNew{5} = '\theta_{8}';

PfAllValueNew(6) = PfAllValue(8);
PfStateStrNew{6} = 'v_{dc8}';

PfAllValueNew(7) = PfAllValue(9);
PfStateStrNew{7} = 'Others';

pie(PfAllValueNew,PfStateStrNew);
if EnableSave
    print(gcf,'Simulation/Figure/InteractionModeParticipation.png','-dpng','-r600');
end