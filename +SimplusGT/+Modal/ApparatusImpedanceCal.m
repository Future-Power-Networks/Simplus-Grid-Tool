% This function calculates the impedance
%
% Author(s): Yue Zhu
% Modified by: Yitong Li, Qipeng Zheng
%
% Notes:
% The unit of lambda is in rad/s

function ZmVal = ApparatusImpedanceCal(GmDSS_Cell, Lambda, ApparatusType)
   
    if ApparatusType <= 89  % AC apparatus
        GmTf.dd=evalfr(GmDSS_Cell(1,1),Lambda);
        GmTf.dq=evalfr(GmDSS_Cell(1,2),Lambda);
        GmTf.qd=evalfr(GmDSS_Cell(2,1),Lambda);
        GmTf.qq=evalfr(GmDSS_Cell(2,2),Lambda);

        Gm = [GmTf.dd, GmTf.dq ; GmTf.qd, GmTf.qq];
        Zm = inv(Gm);

        % Calculate the value of impedance
        ZmVal = Zm;

    elseif ApparatusType >= 1010 && ApparatusType <= 1089 % DC apparatuses  

        GmTf.dd=evalfr(GmDSS_Cell(1,1),Lambda);
        Gm = [GmTf.dd];
        Zm = inv(Gm);
        ZmVal.dd = Zm(1,1);
        ZmVal.dq=[];
        ZmVal.qd=[];
        ZmVal.qq=[];
        ZmVal.d_dc=[];
        ZmVal.q_dc=[];
        ZmVal.dc_d=[];
        ZmVal.dc_q=[];
        ZmVal.dc_dc=[];
        
        % ZmVal.AC_dd = [];
        % ZmVal.AC_dq = [];
        % ZmVal.AC_qd = [];
        % ZmVal.AC_qq = [];
        % 
        % ZmVal.DC = [];

    elseif ApparatusType >= 2000 && ApparatusType <= 2009 % IC apparatuses 

        GmTf.dd=evalfr(GmDSS_Cell(1,1),Lambda);
        GmTf.dq=evalfr(GmDSS_Cell(1,2),Lambda);
        GmTf.qd=evalfr(GmDSS_Cell(2,1),Lambda);
        GmTf.qq=evalfr(GmDSS_Cell(2,2),Lambda);
        GmTf.d_dc=evalfr(GmDSS_Cell(1,3),Lambda);
        GmTf.q_dc=evalfr(GmDSS_Cell(2,3),Lambda);
        GmTf.dc_d=evalfr(GmDSS_Cell(3,1),Lambda);
        GmTf.dc_q=evalfr(GmDSS_Cell(3,2),Lambda);
        GmTf.dc_dc=evalfr(GmDSS_Cell(3,3),Lambda);
        Gm_IC = [GmTf.dd, GmTf.dq ,GmTf.d_dc; GmTf.qd, GmTf.qq, GmTf.q_dc;GmTf.dc_d, GmTf.dc_q, GmTf.dc_dc];
        Gm_AC = [GmTf.dd, GmTf.dq ; GmTf.qd, GmTf.qq];
        Gm_DC = [GmTf.dc_dc];
        Zm_IC = inv(Gm_IC);
        Zm_AC = inv(Gm_AC);
        Zm_DC = inv(Gm_DC);

        % Calculate the value of impedance
        ZmVal.dd = Zm_IC(1,1);
        ZmVal.dq = Zm_IC(1,2);
        ZmVal.qd = Zm_IC(2,1);
        ZmVal.qq = Zm_IC(2,2);
        ZmVal.d_dc = Zm_IC(1,3);
        ZmVal.q_dc = Zm_IC(2,3);
        ZmVal.dc_d = Zm_IC(3,1);
        ZmVal.dc_q = Zm_IC(3,2);
        ZmVal.dc_dc = Zm_IC(3,3);

        % ZmVal.AC_dd = Zm_AC(1,1);
        % ZmVal.AC_dq = Zm_AC(1,2);
        % ZmVal.AC_qd = Zm_AC(2,1);
        % ZmVal.AC_qq = Zm_AC(2,2);
        % 
        % ZmVal.DC = Zm_DC(1,1);
    
    else
        ZmVal.dd=[];
        ZmVal.dq=[];
        ZmVal.qd=[];
        ZmVal.qq=[];
        ZmVal.d_dc=[];
        ZmVal.q_dc=[];
        ZmVal.dc_d=[];
        ZmVal.dc_q=[];
        ZmVal.dc_dc=[];

        % ZmVal.AC_dd = [];
        % ZmVal.AC_dq = [];
        % ZmVal.AC_qd = [];
        % ZmVal.AC_qq = [];
        % 
        % ZmVal.DC = [];
    end
end