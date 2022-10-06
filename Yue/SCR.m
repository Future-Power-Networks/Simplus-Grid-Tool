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
%YADss=GmDss(PortI, PortV);
%YAss=SimplusGT.dss2ss(YADss); %YA
%YAtf=tf(YAss);
%% conventinal SCR, passive components only
Zn2=Zn0(3:4,3:4);
T = [1,1i;1,-1i];% Transform matrix from transfer function to complex vector
Zn2_dqc=T*Zn2*T^(-1);
SCR_1=1/abs(Zn2_dqc(1,1));
SCR_2=1/norm(Zn2_dqc);
fprintf('Conventional SCR at bus-2 is %f. \n',SCR_1);
fprintf('Conventional SCR at bus-2 is calculated from dq-complex norm %f. \n',SCR_2);
%disp('The dq impedance of the system at bus-2 (only passive components) is:');
%disp(Zn2);
disp('The complex dq impedance of the system at bus-2 (only passive components) is:');
disp(Zn2_dqc);
fprintf('===========================================================================\n');
%% SCR, if considering all other components' impedance.
YA0ex2=YA0; YA0ex2(3:4, 3:4) = 0;
Zn0_= (YA0ex2+YN0)^(-1);
Zn02=Zn0_(3:4,3:4);
Zn02_dqc=T*Zn02*T^(-1);
SCR_3= 1/norm(Zn02);
SCR_4= 1/norm(Zn02_dqc);
fprintf('New SCR (including all other apparatus) at bus-2 is calculated from dq-complex norm %f. \n\n',SCR_4);
%disp('The dq impedance of the system at bus-2 (including all other apparatus) is:');
%disp(Zn02);
disp('The complex dq impedance of the system at bus-2 (including all other apparatus) is:');
disp(Zn02_dqc);
fprintf('===========================================================================\n');

SCR_5=1/sqrt(Zn02(1,2)^2+Zn02(2,2)^2);
fprintf('SCR5 is %f. \n\n',SCR_5);


SCR_6=1/sqrt(Zn02(1,1)^2+Zn02(2,1)^2);
fprintf('SCR6 is %f. \n\n',SCR_6);


