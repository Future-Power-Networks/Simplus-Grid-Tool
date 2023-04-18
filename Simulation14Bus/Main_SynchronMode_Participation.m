% Author(s): Yitong Li

clc
close all

EnableSave = 1;

FigN = 0;

FigN = FigN+1;
figure(FigN);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.2 0.4]);
[PfValue,PfStr,Mode] = ModalAnalysisDss(ObjGsysSs,96,6,FigN,0);
% -0.341043723373696 + 50.012489163568944i
PfValueNew = zeros(1,3);
PfStrNew = {'id7-9','iq7-9','Others'};
for i = 1:length(PfValue)
    if strcmp(PfStr(i),PfStrNew{1})
        PfValueNew(1) = PfValueNew(1)+PfValue(i);
    elseif strcmp(PfStr(i),PfStrNew{2})
        PfValueNew(2) = PfValueNew(2)+PfValue(i);
    elseif strcmp(PfStr(i),PfStrNew{3})
        PfValueNew(3) = PfValueNew(3)+PfValue(i);
    else
        error;
    end
end
for i = 1:length(PfStrNew)
    PfStrNew{i} = strrep(PfStrNew{i},'id7-9','i_{d7-9}');
    PfStrNew{i} = strrep(PfStrNew{i},'iq7-9','i_{q7-9}');
end
pie(PfValueNew,PfStrNew);
Mode1 = Mode;
if EnableSave
    print(gcf,'Simulation/Figure/50HzModeSmallDamping.png','-dpng','-r600');
end

FigN = FigN+1;
figure(FigN);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.2 0.4]);
[PfValue,PfStr,Mode] = ModalAnalysisDss(ObjGsysSs,132,2,FigN,0);
PfValueNew = zeros(1,3);
PfStrNew = {'id6-13','iq6-13','Others'};
for i = 1:length(PfValue)
    if strcmp(PfStr(i),PfStrNew{1})
        PfValueNew(1) = PfValueNew(1)+PfValue(i);
    elseif strcmp(PfStr(i),PfStrNew{2})
        PfValueNew(2) = PfValueNew(2)+PfValue(i);
    elseif strcmp(PfStr(i),PfStrNew{3})
        PfValueNew(3) = PfValueNew(3)+PfValue(i);
    else
        error;
    end
end
for i = 1:length(PfStrNew)
    PfStrNew{i} = strrep(PfStrNew{i},'id6-13','i_{d6-13}');
    PfStrNew{i} = strrep(PfStrNew{i},'iq6-13','i_{q6-13}');
end
pie(PfValueNew,PfStrNew);
Mode2 = Mode;
if EnableSave
    print(gcf,'Simulation/Figure/50HzModeLargeDamping.png','-dpng','-r600');
end