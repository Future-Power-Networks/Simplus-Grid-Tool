function [AxisSel, ApparatusSel12, ModeSel, ApparatusSel3,StateSel_DSS, ModeSel_DSS]...
    = ExcelRead(filename, N_Bus, ApparatusType, GminSS)

%% State-PF sheet reading
StateConfig = xlsread(filename,1);
StateSel_DSS = 0;
ModeSel_DSS = 0;
SelIndex1 = 1;
SelIndex2 = 1;

for Count = 1: size(StateConfig,1)
    if StateConfig(Count,1) ==1
        StateSel_DSS(SelIndex1) = Count;
        SelIndex1 = SelIndex1+1;
    end
    if StateConfig(Count,5) ==1
        ModeSel_DSS(SelIndex2) = Count;
        SelIndex2 = SelIndex2+1;
    end
end

%% Impedance-PF sheet reading
%read for impedance-PF.
Config=xlsread(filename,2);
SelIndex = 1;
AxisSel=0;
for ax=1:4
    if Config(ax,4)==1
        AxisSel(SelIndex) = ax;
        SelIndex = SelIndex +1;
    end
end
% if AxisSel == 0
%     error('None axis is selected for bodeplot.');
% end
SelIndex = 1;
ApparatusSel12 = 0;
ApparatusIndex = 1;
for k = 1:N_Bus
    if ApparatusType{k} <= 89 %apparatuses)
        if Config(ApparatusIndex,1) == 1 %the apparatus is selected
            ApparatusSel12(SelIndex) = k;
            SelIndex = SelIndex +1;
        else % not selected.
        end
        ApparatusIndex = ApparatusIndex +1;
    else % floating bus, infinite bus...
    end        
end
% if ApparatusSel12 == 0
%     error('None Apparatus is selected for layer 1&2.');
% end
ModeNum = length(GminSS.A);
SelIndex = 1;
ModeSel =0;
for i=1:ModeNum
    if  Config(i,8)==1
        ModeSel(SelIndex) = i;
        SelIndex = SelIndex+1;
    end
end
% if ModeSel == 0
%     error('None Mode is selected.');
% end
SelIndex = 1;
ApparatusSel3 = 0;
ApparatusIndex = 1;
for k = 1:N_Bus
    if ApparatusType{k} <= 89 %apparatuses)
        if Config(ApparatusIndex,11) == 1 %the apparatus is selected
            ApparatusSel3(SelIndex) = k;
            SelIndex = SelIndex +1;
        else % not selected.
        end
        ApparatusIndex = ApparatusIndex +1;
    else % floating bus, infinite bus...
    end        
end

end