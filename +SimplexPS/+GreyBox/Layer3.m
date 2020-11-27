clear Layer3Result
DeviceSelL3 = 1;
ModeSelL3 = 62;
FreqSelL3=imag(GbMode(ModeSelL3)); 

ZmValOrig = SimplexPS.GreyBox.DeviceImpedanceCal(GmDSS_Cell{DeviceSelL3}, FreqSelL3, DeviceSelL3);
ParamName = fieldnames(Para{DeviceSelL3});
ParamNum = length(ParamName);
ParamValOrig = Para{DeviceSelL3};

%perturb the parameters one by one.
for k=1:ParamNum
    ParaNew = Para;
    ParaSel = getfield(Para{DeviceSelL3},ParamName{k});
    ParaPerturb = ParaSel * 1.01; % 1% perturabation
    ParaNew = setfield(ParaNew{DeviceSelL3}, ParamName{k}, ParaPerturb);
    [~,GmDSS_Cell_New,~,~,~,~,~,~] = ...
        SimplexPS.Toolbox.DeviceModelCreate('Type', DeviceType{DeviceSelL3} ,...
        'Flow',PowerFlow{DeviceSelL3},'Para',ParaNew,'Ts',Ts);
    ZmValNew = SimplexPS.GreyBox.DeviceImpedanceCal(GmDSS_Cell_New, FreqSelL3, DeviceSelL3);
    
    Layer3Result(k).ParaName = ParamName{k};
    Layer3Result(k).DeltaZ.dd = ZmValNew.dd - ZmValOrig.dd;
    Layer3Result(k).DeltaZ.dq = ZmValNew.dq - ZmValOrig.dq;
    Layer3Result(k).DeltaZ.qd = ZmValNew.qd - ZmValOrig.qd;
    Layer3Result(k).DeltaZ.qq = ZmValNew.qq - ZmValOrig.qq;
    
end
