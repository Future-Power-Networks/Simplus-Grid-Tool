function [AxisSel, DeviceSel, ModeSel] = ExcelReadLayer12(filename, N_Bus, DeviceType, GminSS)
Layer12Sel=xlsread(filename,2);
if size(Layer12Sel,2) < 9 % check the size.
    error('Axis, or Device, or Mode is not selected');
end
SelIndex = 1;
AxisSel=0;
for ax=1:4
    if Layer12Sel(ax,4)==1
        AxisSel(SelIndex) = ax;
        SelIndex = SelIndex +1;
    end
end
if AxisSel == 0
    error('None axis is selected.');
end
SelIndex = 1;
DeviceSel = 0;
DeviceIndex = 1;
for k = 1:N_Bus
    if DeviceType{k} <= 89 %devices)
        if Layer12Sel(DeviceIndex,1) == 1 %the device is selected
            DeviceSel(SelIndex) = k;
            SelIndex = SelIndex +1;
        else % not selected.
        end
        DeviceIndex = DeviceIndex +1;
    else % floating bus, infinite bus...
    end        
end
if DeviceSel == 0
    error('None Device is selected.');
end
ModeNum = length(GminSS.A);
SelIndex = 1;
ModeSel =0;
for i=1:ModeNum
    if Layer12Sel(i,9)==1
        ModeSel(SelIndex) = i;
        SelIndex = SelIndex+1;
    end
end
if ModeSel == 0
    error('None Mode is selected.');
end
end