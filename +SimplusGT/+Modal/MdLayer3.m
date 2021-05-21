%Author: Yue Zhu
%Introduction: in this file, parameters of an apparatus Ak is pertrubed one by
%one. The perturbation amplitude Delta_rho = 1e-5*(1+abs(rho)).
% At this toolbox's version, the apparatus's parameters will not affect the 
% results of power flow, therefore, there is no need to calculate power flow again.
% Eventhough the change of rho will affect the equibilium of Ak, it will
% not affect the equilibrium point at other points. Therefore, to calculate
% the new impedance, we only need to update the parameter, and then call the function
% SimplusGT.Toolbox.ApparatusModelCreate, and take the output GmDSS_Cell.
function Layer3Result = MdLayer3(Residue,ZmVal,FreqSel,ApparatusType,...
                ApparatusSelL3All,Para,PowerFlow,Ts,ApparatusBus,ListBus)

ApparatusSelNum=length(ApparatusSelL3All);

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
%         [~,GmDSS_Cell_New,~,~,~,~,~,~] = ...
%             SimplusGT.Toolbox.ApparatusModelCreate_old('Type', ApparatusType{ApparatusSelL3} ,...
%             'Flow',PowerFlow{ApparatusSelL3},'Para',ParaNew,'Ts',Ts);        
        [~,GmDSS_Cell_New,~,~,~,~,~,~,~] ...
        = SimplusGT.Toolbox.ApparatusModelCreate(ApparatusBus{ApparatusSelL3},ApparatusType{ApparatusSelL3},...
                            PowerFlow{ApparatusSelL3},ParaNew,Ts,ListBus);
     
        ZmValNew = SimplusGT.Modal.ApparatusImpedanceCal(GmDSS_Cell_New, FreqSel, ApparatusSelL3);
        
        Layer3Result(ApparatusCount).Apparatus={['Apparatus',num2str(ApparatusSelL3)]};
        Layer3Result(ApparatusCount).Result(k).ParaName = {ParamName{k}};
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
    end
end
end