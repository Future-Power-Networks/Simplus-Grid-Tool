function Layer3Result = GbLayer3(Residue,ZmVal,FreqSel,DeviceType,DeviceSelL3All,Para,PowerFlow,Ts)

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
        ParaSel = getfield(Para{DeviceSelL3},ParamName{k});
        ParaPerturb = ParaSel * 1.01; % 1% perturabation
        ParaNew = setfield(ParaNew{DeviceSelL3}, ParamName{k}, ParaPerturb);
        [~,GmDSS_Cell_New,~,~,~,~,~,~] = ...
            SimplexPS.Toolbox.DeviceModelCreate('Type', DeviceType{DeviceSelL3} ,...
            'Flow',PowerFlow{DeviceSelL3},'Para',ParaNew,'Ts',Ts);
        ZmValNew = SimplexPS.GreyBox.DeviceImpedanceCal(GmDSS_Cell_New, FreqSel, DeviceSelL3);

        Layer3Result{DeviceCount}(k).ParaName = ParamName{k};
        Layer3Result{DeviceCount}(k).DeltaZ.dd = ZmValNew.dd - ZmValOrig.dd;
        Layer3Result{DeviceCount}(k).DeltaZ.dq = ZmValNew.dq - ZmValOrig.dq;
        Layer3Result{DeviceCount}(k).DeltaZ.qd = ZmValNew.qd - ZmValOrig.qd;
        Layer3Result{DeviceCount}(k).DeltaZ.qq = ZmValNew.qq - ZmValOrig.qq;
        
         Layer3Result{DeviceCount}(k).DLambda = -1*(...
             Layer3Result{DeviceCount}(k).DeltaZ.dd * Residue_.dd...
            + Layer3Result{DeviceCount}(k).DeltaZ.dq * Residue_.qd ...
            + Layer3Result{DeviceCount}(k).DeltaZ.qd * Residue_.dq ...
            + Layer3Result{DeviceCount}(k).DeltaZ.qq * Residue_.qq);
    end
end
end