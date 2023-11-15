% This function calculates the impedance
%
% Author(s): Yue Zhu
% Modified by: Yitong Li, Qipeng Zheng
%
% Notes:
% The unit of lambda is in rad/s

function ZmVal = ApparatusImpedanceCal(GmDSS_Cell, Lambda, ApparatusType)
   
    if ApparatusType <= 89  % Ac apparatus
        GmTf = evalfr(GmDSS_Cell(1:2,1:2),Lambda);
        ZmVal = inv(GmTf);
    elseif ApparatusType >= 1000 && ApparatusType <= 1089 % Dc apparatus
        GmTf = evalfr(GmDSS_Cell(1:1,1:1),Lambda);
        ZmVal = inv(GmTf);
    elseif ApparatusType >= 2000 && ApparatusType <= 2009 % Interlink apparatus
        GmTf = evalfr(GmDSS_Cell(1:3,1:3),Lambda);
        ZmVal = inv(GmTf);
    else
        ZmVal = [];
    end

end