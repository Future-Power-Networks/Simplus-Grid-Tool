ZmDSS_cell = cell (N_Bus);
ZmDSS = dss([],[],[],[],[]);
for k = 1:N_Bus
    YmObj_ = GmObj_Cell{k};
    [~,YmDSS_] = YmObj_.GetDSS(YmObj_);
    if isempty(YmDSS_.A)
        ZmDSS_cell{k} = YmDSS_;
    else
        YmDSS_ = YmDSS_(1:2,1:2);
        ZmDSS_cell{k} = inv(YmDSS_);
    end
    ZmDSS = SimplusGT.DssAppend(ZmDSS,ZmDSS_cell{1});
end
ZsysDSS = feedback(ZmDSS,YbusDSS);
[za,zb] = tfdata(ZsysDSS(15,15));
%minreal(ZsysDSS)
[~,ZbusDSS] = ZbusObj.GetDSS(ZbusObj);
% Ym8Obj = GmObj_Cell{8};
% [~,Ym8DSS] = Ym8Obj.GetDSS(Ym8Obj);
% Ym8_dss = Ym8DSS(1:2,1:2);
% Zm8DSS = inv(Ym8_dss);
% %Zm8tfdd = tf(c);
% [ znum , zden ] = tfdata( Zm8DSS(1,1) );
% [zr,zp,~] = residue(znum{1},zden{1});
% zr = zr/2/pi;
%tf_Ym8=zeros(2);
% tf_Ym8(1,1) = tf(Ym8DSS(1,1));
% tf_Ym8(1,2) = tf(Ym8DSS(1,2));
% tf_Ym8(2,1) = tf(Ym8DSS(2,1));
% tf_Ym8(2,2) = tf(Ym8DSS(2,2));
% Ym8SS = tf2ss(

%Ym8DSS.B = Ym8DSS.B(:,1:2);
%Ym8DSS.D= Ym8DSS.D(1:2, 1:2);
%Ym8DSS.C= Ym8DSS.C(1:2,:);

%Ym8DSS.B = Bx;

%inv(Ym8DSS)h