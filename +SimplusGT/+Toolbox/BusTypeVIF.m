% Find the voltage, current, floating buses
%
% Author(s): Jialu Yuan, Yitong Li
%
% Notes:
% Voltage bus (vbus): GFM inverter, SG, inf bus, etc.
% Current bus (ibus): GFL inverter, etc.
% Floating bus (fbus): floating node, passive node, etc.
% 
% Essentially, the floating bus can also be regarded as ibus with zero current.

function [vbus,ibus,fbus] = BusTypeVIF(ApparatusType)

l = 1;
m = 1;
n = 1;

vbus = [];
ibus = [];
fbus = [];
for i = 1:length(ApparatusType)
    if ((ApparatusType{i} >= 0) && (ApparatusType{i} <= 9)) || (ApparatusType{i} == 90) || ((ApparatusType{i} >= 20) && (ApparatusType{i} <= 40)) || (ApparatusType{i} == 50)
        vbus(l) = i;    
        l = l+1;
    elseif ((ApparatusType{i} >= 10) && (ApparatusType{i} <= 19)) || (ApparatusType{i} == 41) || (ApparatusType{i} == 51)
        ibus(m) = i;   
        m = m+1;
    elseif (ApparatusType{i} == 100)
        fbus(n) = i;
        n = n+1;
    else
        error('Error');
    end
end

end