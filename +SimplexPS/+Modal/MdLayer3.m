%Author: Yue Zhu
%Introduction: in this file, parameters of an apparatus Ak is pertrubed one by
%one. The perturbation amplitude Delta_rho = 1e-5*(1+abs(rho)).
% At this toolbox's version, the apparatus's parameters will not affect the 
% results of power flow, therefore, there is no need to calculate power flow again.
% Eventhough the change of rho will affect the equibilium of Ak, it will
% not affect the equilibrium point at other points. Therefore, to calculate
% the new impedance, we only need to update the parameter, and then call the function
% SimplexPS.Toolbox.DeviceModelCreate, and take the output GmDSS_Cell.

function Layer3Result = MdLayer3(Residue,ZmVal,FreqSel,DeviceType,DeviceSelL3All,Para,PowerFlow,Ts)

DeviceSelNum=length(DeviceSelL3All);

for DeviceCount = 1:DeviceSelNum
    DeviceSelL3 = DeviceSelL3All(DeviceCount);
    ZmValOrig = ZmVal(DeviceSelL3);
    ParamName = fieldnames(Para{DeviceSelL3});
    ParamNum = length(ParamName);
    Residue_ = Residue(DeviceSelL3);
    %perturb the parameters one by one.
    for k=1:ParamNum
        ParaNew = Para;
        ParaSel = getfield(Para{DeviceSelL3},ParamName{k}); % extract the parameter value
        delta_para = 1e-5*(1+abs(ParaSel));
        ParaPerturb = ParaSel + delta_para ; % add perturabation
        ParaNew = setfield(ParaNew{DeviceSelL3}, ParamName{k}, ParaPerturb); % update the parameter       
        [~,GmDSS_Cell_New,~,~,~,~,~,~] = ...
            SimplexPS.Toolbox.DeviceModelCreate('Type', DeviceType{DeviceSelL3} ,...
            'Flow',PowerFlow{DeviceSelL3},'Para',ParaNew,'Ts',Ts);
        ZmValNew = SimplexPS.Modal.DeviceImpedanceCal(GmDSS_Cell_New, FreqSel, DeviceSelL3);
        
        Layer3Result(DeviceCount).Device={['Device',num2str(DeviceSelL3)]};
        %Layer3Result(DeviceCount).Result(k)={['Device',num2str(DeviceSelL3)]};
        Layer3Result(DeviceCount).Result(k).ParaName = {ParamName{k}};
        Layer3Result(DeviceCount).Result(k).DeltaZ.dd = (ZmValNew.dd - ZmValOrig.dd)/(delta_para);
        Layer3Result(DeviceCount).Result(k).DeltaZ.dq = (ZmValNew.dq - ZmValOrig.dq)/(delta_para);
        Layer3Result(DeviceCount).Result(k).DeltaZ.qd = (ZmValNew.qd - ZmValOrig.qd)/(delta_para);
        Layer3Result(DeviceCount).Result(k).DeltaZ.qq = (ZmValNew.qq - ZmValOrig.qq)/(delta_para);
        
        Layer3Result(DeviceCount).Result(k).Residue.dd = Residue_.dd;
        Layer3Result(DeviceCount).Result(k).Residue.dq = Residue_.dq;
        Layer3Result(DeviceCount).Result(k).Residue.qd = Residue_.qd;
        Layer3Result(DeviceCount).Result(k).Residue.qq = Residue_.qq;
        
         Layer3Result(DeviceCount).Result(k).DLambda_rad = -1*(...
             Layer3Result(DeviceCount).Result(k).DeltaZ.dd * Residue_.dd...
            + Layer3Result(DeviceCount).Result(k).DeltaZ.dq * Residue_.qd ...
            + Layer3Result(DeviceCount).Result(k).DeltaZ.qd * Residue_.dq ...
            + Layer3Result(DeviceCount).Result(k).DeltaZ.qq * Residue_.qq);
        DLambda_Hz=Layer3Result(DeviceCount).Result(k).DLambda_rad/(2*pi);
        Layer3Result(DeviceCount).Result(k).DLambdaRho_Hz=DLambda_Hz;
        Layer3Result(DeviceCount).Result(k).DLambdaRho_pu_Hz=DLambda_Hz*ParaSel;
    end
end
end