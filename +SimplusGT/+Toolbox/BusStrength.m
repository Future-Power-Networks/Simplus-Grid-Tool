% Calculate the Thevenin impedance seen at each bus of a power system. 
%
% Author(s): Jialu Yuan, Yitong Li

function [Ydiag,Ybus] = BusStrength(ApparatusType,ListLine)

[vbus,ibus,fbus] = SimplusGT.Toolbox.BusTypeVIF(ApparatusType);

% Calculate Ybus
Ybus = SimplusGT.PowerFlow.YbusCalc(ListLine); 
 
% Sort Ybus
Ybus_1 = Ybus([vbus,ibus,fbus],:);
YbusSort = Ybus_1(:,[vbus,ibus,fbus]);

% Calculate Hbus
HbusSort = SimplusGT.Strength.Y2H(YbusSort,length(vbus));

% Sort Hbus
Hbus_1 ([vbus,ibus,fbus],:) = HbusSort(:,:);
Hbus(:,[vbus,ibus,fbus]) = Hbus_1(:,:);

% Calculate bus impedance
for i = 1:length(Hbus)
    if (find(vbus,i))
        Ydiag(i) = Hbus(i,i);
    else
        Ydiag(i) = 1/Hbus(i,i);
    end
end

end

