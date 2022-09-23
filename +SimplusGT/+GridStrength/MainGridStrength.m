% Author(s): Yitong

% Get Zm
ObjYmDss = SimplusGT.ObjTruncate(ObjGmDss,PortI,PortV);
[~,~,LengthYmDss] = SimplusGT.ObjSize(ObjYmDss);
ObjZmDss = SimplusGT.ObjSwitchInOut(ObjYmDss,LengthYmDss);
clear('LengthYmDss')

% Connect Zm and Ybus
ObjZsysDss = SimplusGT.ObjFeedback(ObjZmDss,ObjYbusDss);

[~,ZsysDss] = ObjZsysDss.GetDSS(ObjZsysDss);
ZsysSs = minreal(ZsysDss);

stop