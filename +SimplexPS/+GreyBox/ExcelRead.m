function [AxisSel, DeviceSel12, ModeSel, DeviceSel3] = ExcelRead(filename, N_Bus, DeviceType, GminSS)
%filename='GreyBoxConfig.xlsx';
Config=xlsread(filename,2);


%%
%read for impedance-PF.
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
DeviceSel12 = 0;
DeviceIndex = 1;
for k = 1:N_Bus
    if DeviceType{k} <= 89 %devices)
        if Config(DeviceIndex,1) == 1 %the device is selected
            DeviceSel12(SelIndex) = k;
            SelIndex = SelIndex +1;
        else % not selected.
        end
        DeviceIndex = DeviceIndex +1;
    else % floating bus, infinite bus...
    end        
end
% if DeviceSel12 == 0
%     error('None Device is selected for layer 1&2.');
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
DeviceSel3 = 0;
DeviceIndex = 1;
for k = 1:N_Bus
    if DeviceType{k} <= 89 %devices)
        if Config(DeviceIndex,11) == 1 %the device is selected
            DeviceSel3(SelIndex) = k;
            SelIndex = SelIndex +1;
        else % not selected.
        end
        DeviceIndex = DeviceIndex +1;
    else % floating bus, infinite bus...
    end        
end

end