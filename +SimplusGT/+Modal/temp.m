 mode_rad = pole_sys(186)*2*pi;
% Gm = GmDSS_Cell{3};
% Gm = Gm(1:2,1:2);
% Zm = inv(Gm);
% ZmValOrig.dd=evalfr(Zm(1,1),mode_rad);
% ZmValOrig.dq=evalfr(Zm(1,2),mode_rad);
% ZmValOrig.qd=evalfr(Zm(2,1),mode_rad);
% ZmValOrig.qq=evalfr(Zm(2,2),mode_rad);
% 
% ParamName = fieldnames(Para{3});
% ParaNew = Para;
% ParaSel = getfield(Para{3},ParamName{1}); % extract the parameter
% delta_para = 0;%1e-5*(1+abs(ParaSel));
% ParaPerturb = ParaSel + delta_para ; % add perturabation
% ParaNew = setfield(ParaNew{3}, ParamName{1}, ParaPerturb); % update the parameter  
% [~,GmDSS_Cell_New,~,~,~,~,~,~,~] ...
%         = SimplusGT.Toolbox.ApparatusModelCreate(ApparatusBus{3},ApparatusType{3},...
%                             PowerFlow{3},Para{3},Ts,ListBus);
% 
% Gm = GmDSS_Cell_New;
% Gm = Gm(1:2,1:2);
% Zm = inv(Gm);
% ZmValNew.dd=evalfr(Zm(1,1),mode_rad);
% ZmValNew.dq=evalfr(Zm(1,2),mode_rad);
% ZmValNew.qd=evalfr(Zm(2,1),mode_rad);
% ZmValNew.qq=evalfr(Zm(2,2),mode_rad);


GmDSS_Cell_3=GmDSS_Cell{3};
GmValOrig.dd=evalfr(GmDSS_Cell_3(1,1),mode_rad)

[~,GmDSS_Cell_3x,~,~,~,~,~,~,~] ...
        = SimplusGT.Toolbox.ApparatusModelCreate(ApparatusBus{3},ApparatusType{3},...
                            ApparatusPowerFlow{3},Para{3},Ts,ListBus);
                        
GmValOrig.dd=evalfr(GmDSS_Cell_3x(1,1),mode_rad)



% function [GmObj,GmDSS,ApparatusPara,ApparatusEqui,DiscreDampingResistor,OtherInputs,StateStr,InputStr,OutputStr] ...
%         = ApparatusModelCreate(ApparatusBus,Type,PowerFlow,Para,Ts,ListBus) 
%     
% [GmObj_Cell{i},GmDSS_Cell{i},ApparatusPara{i},ApparatusEqui{i},ApparatusDiscreDamping{i},OtherInputs{i},ApparatusStateStr{i},ApparatusInputStr{i},ApparatusOutputStr{i}] = ...
%         SimplusGT.Toolbox.ApparatusModelCreate(ApparatusBus{i},ApparatusType{i},ApparatusPowerFlow{i},Para{i},Ts,ListBus);
    