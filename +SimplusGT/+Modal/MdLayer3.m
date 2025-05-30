% Layer 3 modal analysis.
%
% Author(s): Yue Zhu
% Modified by: Yitong Li, Qipeng Zheng
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
    ZmValOrig = ZmVal{ApparatusSelL3};  
    ParamName = fieldnames(Para{ApparatusSelL3});
    ParamNum = length(ParamName);
    Residue_ = Residue{ApparatusSelL3};
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
     
        ZmValNew = SimplusGT.Modal.ApparatusImpedanceCal(GmDSS_Cell_New, Mode_rad, ApparatusType{ApparatusSelL3});
        
        Layer3Result(ApparatusCount).Apparatus={['Apparatus',num2str(ApparatusSelL3)]};
        Layer3Result(ApparatusCount).Result(k).ParaName = ParamName(k);

        Layer3Result(ApparatusCount).Result(k).DeltaZ = (ZmValNew-ZmValOrig)/(delta_para);
        % conj(sum(dot(A,B'))) = A(1,1)*B(1,1) + A(1,2)*B(2,1) + A(2,1)*B(1,2) + A(2,2)*B(2,2)
        Layer3Result(ApparatusCount).Result(k).DLambda_rad = -1*conj(sum(dot(Residue_,Layer3Result(ApparatusCount).Result(k).DeltaZ')));
        DLambda_Hz=Layer3Result(ApparatusCount).Result(k).DLambda_rad/(2*pi);
        Layer3Result(ApparatusCount).Result(k).DLambdaRho_Hz=DLambda_Hz;
        Layer3Result(ApparatusCount).Result(k).DLambdaRho_pu_Hz=DLambda_Hz*ParaSel;

    end
end
end