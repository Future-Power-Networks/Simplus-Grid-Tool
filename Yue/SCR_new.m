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
kp=find(OrderOld2New==Poi); % Number of the interested POI in the hybrid matrix
%% SCR1:  conventional SCR: considering all sources as voltage source 
if ApparatusSourceType(kp)==1 % voltage source
    YnetP = HN(kp*2-1: kp*2, kp*2-1: kp*2);
    ZnetP = inv(YnetP);
elseif ApparatusSourceType(kp)==2 % current source
    ZnetP = HN(kp*2-1: kp*2, kp*2-1: kp*2);
    %YnetP = inv(ZnetP);
end
SCR1 = 1/abs_z(ZnetP);
fprintf('Conventional SCR at bus-%d is %f. \n',[Poi,SCR1]);
%% SCR2: Equivalent Circuit-based Short Circuit Ratio: ESCR
% GFL ONLY
POI_IF=[0,0;0,0]; % interaction factor
for i=1:NumBus
     if ApparatusSourceType(i)==1 % voltage source, GFM
         POI_IF=POI_IF + HN(kp*2-1:kp*2 , i*2-1:i*2);
     elseif ApparatusSourceType(i)==2
         Znet_i = HN(i*2-1:i*2 , i*2-1:i*2); 
         POI_IF=POI_IF + HN(kp*2-1:kp*2 , i*2-1:i*2) / Znet_i ; 
     end
end
ESCR = 1/abs_z(ZnetP * POI_IF);
fprintf('ESCR at bus-%d is %f. \n',[Poi,ESCR]);
%% SCR3: ESCR-NEW: including the dynamics of other apparatus
Ypp=YA0(2*Poi-1:2*Poi, 2*Poi-1:2*Poi);
Zpp=inv(Ypp);
POI_IF_2=[0,0;0,0]; % interaction factor
for i=1:NumBus
     k=OrderOld2New(i); % the original bus number
     if ApparatusSourceType(i)==1 % voltage source, GFM
         YAkk=YA0(2*k-1:2*k, 2*k-1:2*k);
         ZAkk=inv(YAkk);
         POI_IF_2=POI_IF_2 + HN(kp*2-1:kp*2,i*2-1:i*2) + ZAkk/ZnetP;
     elseif ApparatusSourceType(i)==2
         YAkk=YA0(2*k-1:2*k, 2*k-1:2*k);
         Ynetk = inv(HN(kp*2-1:kp*2,i*2-1:i*2));
         Ynet_combine = YAkk+Ynetk;
         Znet_combine = inv(Ynet_combine);
         Znet_i = HN(i*2-1:i*2 , i*2-1:i*2); 
         POI_IF_2=POI_IF_2 + Znet_combine/Znet_i ; 
     end
end
ESCR_new = 1/abs_z(ZnetP * POI_IF_2);
fprintf('ESCR-new at bus-%d is %f. \n',[Poi,ESCR_new]);
%abs_z(X)
function X=abs_z(G)
    X=sqrt( (G(1,1)+G(1,2))^2 + (G(2,1)+G(2,2))^2 )/sqrt(2);
end
