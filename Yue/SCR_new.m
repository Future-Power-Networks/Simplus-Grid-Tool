%% Get all the matrices.
sc=0; % s=0 for d-q frame.
YNdss = YbusDss;%SimplusGT.dss2ss(YbusDss); % YN
ObjZNdss = SimplusGT.ObjDss2Ss(ObjZbusDss);
[~,ZN] = ObjYsysSs.GetSS(ObjZNdss); % Znodal.
ObjYAss=SimplusGT.ObjDss2Ss(ObjYmDss);
[~,YA]= ObjYAss.GetSS(ObjYAss);
Zn0=evalfr(ZN,sc); % get numerical matrix -- nodal impedance matrix
YA0=evalfr(tf(YA),sc);
YN0=evalfr(YbusDss,sc);
% get the hybrid nodal matrix GHm.
[HN, OrderOld2New, ApparatusSourceType] = SimplusGT.GridStrength.HybridMatrix(YN0, ListBus, ApparatusType, NumBus);


Poi=2; % interested POI.
kp=find(OrderOld2New==2); % order of the interested POI.
%% SCR1:  conventional SCR: considering all sources as voltage source 
if ApparatusSourceType(kp)==1 % voltage source
    YnetP = HN(kp*2-1: kp*2, kp*2-1: kp*2);
    ZnetP = inv(YnetP);
elseif ApparatusSourceType(kp)==2 % current source
    ZnetP = HN(kp*2-1: kp*2, kp*2-1: kp*2);
    %YnetP = inv(ZnetP);
end
SCR1 = 1/norm(ZnetP);
fprintf('Conventional SCR at bus-%d is %f. \n',[Poi,SCR1]);
%% SCR2: Equivalent Circuit-based Short Circuit Ratio: ESCR
% GFL ONLY
POI_IF=[0,0;0,0]; % interaction factor
for i=1:NumBus
     if ApparatusSourceType(i)==1 % voltage source, GFM
         POI_IF=POI_IF + HN(kp*2-1:kp*2 , i*2-1:i*2);
     elseif ApparatusSourceType(i)==2
         POI_IF=POI_IF + HN(kp*2-1:kp*2 , i*2-1:i*2) / norm(HN(kp*2-1: kp*2, kp*2-1: kp*2)) ; 
     end
end
ESCR = 1/norm(ZnetP * POI_IF);
fprintf('ESCR at bus-%d is %f. \n',[Poi,ESCR]);
%% SCR4: Weighted short circuit ratio (WSCR)
