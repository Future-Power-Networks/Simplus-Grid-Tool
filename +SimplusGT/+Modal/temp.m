% ZmDSS_cell = cell(N_Bus);
% ZmDSS = dss([],[],[],[],[]);
% for k = 1:N_Bus
%     YmObj_ = GmObj_Cell{k};
%     [~,YmDSS_] = YmObj_.GetDSS(YmObj_);
%     if isempty(YmDSS_.A)    %floating bus
%         ZmDSS_cell{k} = YmDSS_;
%     else
%         YmDSS_ = YmDSS_(1:2,1:2);
%         ZmDSS_cell{k} = inv(YmDSS_);
%         %ZmDSS_cell{k} = SimplusGT.DssSwitchInOut(YmDSS_,2);
%     end
%     ZmDSS = SimplusGT.DssAppend(ZmDSS,ZmDSS_cell{k});
% end

[ZsysObj,ZsysDSS] = SimplusGT.WholeSysZ_cal(GmObj_Cell,YbusObj,N_Apparatus, N_Bus);

%ZsysDSS = feedback(ZmDSS,YbusDSS);
%[za,zb] = tfdata(ZsysDSS(15,15));
Zsys_pole = pole(ZsysDSS)/2/pi;
figure(2);
clf
scatter(real(Zsys_pole)*2*pi,imag(Zsys_pole),'x','LineWidth',1.5); hold on; grid on;
xlabel('Real Part (Hz)');
ylabel('Imaginary Part (Hz)');
title('Zoomed pole map');
axis([-30,10,-3.5,3.5]);

[~,ZmInStr,ZmOutStr] = ZsysObj.GetString(ZsysObj);%get name string
Port_i_in = [];
Port_v_out = [];
for i = 1:N_Apparatus
    [~,in1] = SimplusGT.CellFind(ZmInStr,['i_d',num2str(i)]);    
    [~,out1] = SimplusGT.CellFind(ZmOutStr,['v_d',num2str(i)]);
    if ~isempty(in1) %ac apparatus
        Port_i_in = [Port_i_in,in1,in1+1];
        Port_v_out = [Port_v_out,out1,out1+1];
    else
        error(['Error']);
    end
end


%minreal(ZsysDSS)
%[~,ZbusDSS] = ZbusObj.GetDSS(ZbusObj);
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