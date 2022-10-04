%% Get all the matrices.
YNdss = YbusDss;%SimplusGT.dss2ss(YbusDss); % YN
ObjZNdss = SimplusGT.ObjDss2Ss(ObjZbusDss);
[~,ZN] = ObjYsysSs.GetSS(ObjZNdss); % Znodal.
[~,YA]= ObjYmDss.GetDSS(ObjYmDss);
%YADss=GmDss(PortI, PortV);
%YAss=SimplusGT.dss2ss(YADss); %YA
%YAtf=tf(YAss);
%% conventinal SCR, passive components only
ZN2_dq=[ZN(3,3),ZN(3,4); ZN(4,3),ZN(4,4)];
T = [1,1i;1,-1i];% Transform matrix from transfer function to complex vector
ZN2_dqc=T*ZN2_dq*T^(-1);
ZN2_dqc_0=evalfr(ZN2_dqc,0);
SNR_n2_1=1/abs(ZN2_dqc_0(1,1));
SNR_n2_2=1/norm(ZN2_dqc_0);
fprintf('Conventional SNR at bus-2 is %f. \n',SNR_n2);
fprintf('Conventional SNR at bus-2 is calculated from dq-complex norm %f. \n',SNR_n2_2);

%% SCR, if considering all other components' impedance.
YA_ex2 = YA; YA_ex2.B(:,3:4)=0;  YA_ex2.C(3:4,:)=0; % take out YA2 from YA
YA_ex2_0 = evalfr(tf(YA_ex2),0);
YN_0 = evalfr(YNdss,0);
Znodal_ex2_0 = (YA_ex2_0+YN_0)^(-1);
Znodal2_0=Znodal_ex2_0(3:4,3:4);
Znodal2_dqc_0=T*Znodal2_0*T^(-1);
disp('The impedance of the system at bus-2 (including all other apparatus) is:');
disp(Znodal2_dqc_0);%Znodal2_dqc_0
SNR_n3=1/norm(Znodal2_dqc_0);
fprintf('SNR (including all other apparatus) at bus-2 is calculated from dq-complex norm %f. \n',SNR_n3);

%Ynodal_ex2 = (YAdss_ex2) + (YbusDss); %%%%%% not correct. shouldn't be added directly.


%Znodal_ex2 = SimplusGT.DssSwitchInOut(Ynodal_ex2,14);
%ZN2_dq_ex=[Znodal_ex2(3,3),Znodal_ex2(3,4); Znodal_ex2(4,3),Znodal_ex2(4,4)];
%ZN2_dqc_ex=T*ZN2_dq_ex*T^(-1);
%ZN2_dqc_ex_0=evalfr(ZN2_dqc_ex,0);

