% Layer 3 modal analysis.
%
% Author(s): Yue Zhu
% Modifed by: Yitong Li, Qipeng Zheng
%
%Introduction: in this file, parameters of an apparatus Ak is pertrubed one by
%one. The perturbation amplitude Delta_rho = 1e-5*(1+abs(rho)).
% At this toolbox's version, the apparatus's parameters will not affect the 
% results of power flow, therefore, there is no need to calculate power flow again.
% Eventhough the change of rho will affect the equibilium of Ak, it will
% not affect the equilibrium point at other points. Therefore, to calculate
% the new impedance, we only need to update the parameter, and then call the function
% SimplusGT.Toolbox.ApparatusModelCreate, and take the output GmDSS_Cell.

function Layer3Result = MdLayer3(Residue,ZmVal,Mode_Hz,ApparatusType,...
                ApparatusSelL3All,Para,ApparatusPowerFlow,Ts,ApparatusBus,ListBus)

ApparatusSelNum=length(ApparatusSelL3All);
Mode_rad = Mode_Hz*2*pi;
for ApparatusCount = 1:ApparatusSelNum
    ApparatusSelL3 = ApparatusSelL3All(ApparatusCount);
    ZmValOrig = ZmVal(ApparatusSelL3);  
    ParamName = fieldnames(Para{ApparatusSelL3});
    ParamNum = length(ParamName);
    Residue_ = Residue(ApparatusSelL3);
    %perturb the parameters one by one.
    for k=1:ParamNum
        ParaNew = Para;
        ParaSel = getfield(Para{ApparatusSelL3},ParamName{k}); % extract the parameter
        delta_para = 1e-5*(1+abs(ParaSel));
        ParaPerturb = ParaSel + delta_para ; % add perturabation
        ParaNew = setfield(ParaNew{ApparatusSelL3}, ParamName{k}, ParaPerturb); % update the parameter  
   
        [~,GmDSS_Cell_New,~,~,~,~,~,~,~] ...
        = SimplusGT.Toolbox.ApparatusModelCreate(ApparatusBus{ApparatusSelL3},ApparatusType{ApparatusSelL3},...
                            ApparatusPowerFlow{ApparatusSelL3},ParaNew,Ts,ListBus);
     
       ZmValNew = SimplusGT.Modal.ApparatusImpedanceCal(GmDSS_Cell_New, Mode_rad);
        
        Layer3Result(ApparatusCount).Apparatus={['Apparatus',num2str(ApparatusSelL3)]};
        Layer3Result(ApparatusCount).Result(k).ParaName = ParamName(k);

        if ApparatusType{ApparatusSelL3} <= 89 % Ac apparatus
            Layer3Result(ApparatusCount).Result(k).DeltaZ.dd = (ZmValNew.dd - ZmValOrig.dd)/(delta_para);
            Layer3Result(ApparatusCount).Result(k).DeltaZ.dq = (ZmValNew.dq - ZmValOrig.dq)/(delta_para);
            Layer3Result(ApparatusCount).Result(k).DeltaZ.qd = (ZmValNew.qd - ZmValOrig.qd)/(delta_para);
            Layer3Result(ApparatusCount).Result(k).DeltaZ.qq = (ZmValNew.qq - ZmValOrig.qq)/(delta_para);
            
              Layer3Result(ApparatusCount).Result(k).DLambda_rad = -1*(...
                 Layer3Result(ApparatusCount).Result(k).DeltaZ.dd * Residue_.dd...
                + Layer3Result(ApparatusCount).Result(k).DeltaZ.dq * Residue_.qd ...
                + Layer3Result(ApparatusCount).Result(k).DeltaZ.qd * Residue_.dq ...
                + Layer3Result(ApparatusCount).Result(k).DeltaZ.qq * Residue_.qq);
            DLambda_Hz=Layer3Result(ApparatusCount).Result(k).DLambda_rad/(2*pi);
            Layer3Result(ApparatusCount).Result(k).DLambdaRho_Hz=DLambda_Hz;
            Layer3Result(ApparatusCount).Result(k).DLambdaRho_pu_Hz=DLambda_Hz*ParaSel;

        elseif ApparatusType{ApparatusSelL3} >= 1010 && ApparatusType{ApparatusSelL3} <= 1089 % DC apparatuses
            Layer3Result(ApparatusCount).Result(k).DeltaZ.dd = (ZmValNew.dd - ZmValOrig.dd)/(delta_para);
            Layer3Result(ApparatusCount).Result(k).DLambda_rad = -1*(...
                Layer3Result(ApparatusCount).Result(k).DeltaZ.dd * Residue_.dd);
            DLambda_Hz=Layer3Result(ApparatusCount).Result(k).DLambda_rad/(2*pi);
            Layer3Result(ApparatusCount).Result(k).DLambdaRho_Hz=DLambda_Hz;
            Layer3Result(ApparatusCount).Result(k).DLambdaRho_pu_Hz=DLambda_Hz*ParaSel;

        elseif ApparatusType{ApparatusSelL3} >= 2000 && ApparatusType{ApparatusSelL3} <= 2009 % IC apparatuses
            Layer3Result(ApparatusCount).Result(k).DeltaZ.dd = (ZmValNew.dd - ZmValOrig.dd)/(delta_para);
            Layer3Result(ApparatusCount).Result(k).DeltaZ.dq = (ZmValNew.dq - ZmValOrig.dq)/(delta_para);
            Layer3Result(ApparatusCount).Result(k).DeltaZ.qd = (ZmValNew.qd - ZmValOrig.qd)/(delta_para);
            Layer3Result(ApparatusCount).Result(k).DeltaZ.qq = (ZmValNew.qq - ZmValOrig.qq)/(delta_para);
            Layer3Result(ApparatusCount).Result(k).DeltaZ.d_dc = (ZmValNew.d_dc - ZmValOrig.d_dc)/(delta_para);
            Layer3Result(ApparatusCount).Result(k).DeltaZ.q_dc = (ZmValNew.q_dc - ZmValOrig.q_dc)/(delta_para);
            Layer3Result(ApparatusCount).Result(k).DeltaZ.dc_d = (ZmValNew.dc_d - ZmValOrig.dc_d)/(delta_para);
            Layer3Result(ApparatusCount).Result(k).DeltaZ.dc_q = (ZmValNew.dc_q - ZmValOrig.dc_q)/(delta_para);
            Layer3Result(ApparatusCount).Result(k).DeltaZ.dc_dc = (ZmValNew.dc_dc - ZmValOrig.dc_dc)/(delta_para);

            % Layer3Result(ApparatusCount).Result(k).DeltaZ.AC_dd = (ZmValNew.AC_dd - ZmValOrig.AC_dd)/(delta_para);
            % Layer3Result(ApparatusCount).Result(k).DeltaZ.AC_dq = (ZmValNew.AC_dq - ZmValOrig.AC_dq)/(delta_para);
            % Layer3Result(ApparatusCount).Result(k).DeltaZ.AC_qd = (ZmValNew.AC_qd - ZmValOrig.AC_qd)/(delta_para);
            % Layer3Result(ApparatusCount).Result(k).DeltaZ.AC_qq = (ZmValNew.AC_qq - ZmValOrig.AC_qq)/(delta_para);
            % Layer3Result(ApparatusCount).Result(k).DeltaZ.DC = (ZmValNew.DC - ZmValOrig.DC)/(delta_para);

            Layer3Result(ApparatusCount).Result(k).DLambda_rad_IC = -1*(...
                 Layer3Result(ApparatusCount).Result(k).DeltaZ.dd * Residue_.dd...
                + Layer3Result(ApparatusCount).Result(k).DeltaZ.dq * Residue_.qd ...
                + Layer3Result(ApparatusCount).Result(k).DeltaZ.qd * Residue_.dq ...
                + Layer3Result(ApparatusCount).Result(k).DeltaZ.qq * Residue_.qq...
                + Layer3Result(ApparatusCount).Result(k).DeltaZ.d_dc * Residue_.dc_d ...
                + Layer3Result(ApparatusCount).Result(k).DeltaZ.q_dc * Residue_.dc_q ...
                + Layer3Result(ApparatusCount).Result(k).DeltaZ.dc_d * Residue_.d_dc ...
                + Layer3Result(ApparatusCount).Result(k).DeltaZ.dc_q * Residue_.q_dc ...
                + Layer3Result(ApparatusCount).Result(k).DeltaZ.dc_dc * Residue_.dc_dc);
       
            DLambda_Hz_IC=Layer3Result(ApparatusCount).Result(k).DLambda_rad_IC/(2*pi);
            Layer3Result(ApparatusCount).Result(k).DLambdaRho_Hz_IC=DLambda_Hz_IC;
        end
    end
end
end