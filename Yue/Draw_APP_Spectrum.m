ObjYAss=SimplusGT.ObjDss2Ss(ObjYmDss);
[~,YA]= ObjYAss.GetSS(ObjYAss);
figure(5);
clf;
subplot(2,2,1)
bode(YA(3,3))
hold on;
bode(YA(5,5));
subplot(2,2,2)
bode(YA(3,4))
hold on;
bode(YA(5,6));
subplot(2,2,3)
bode(YA(4,3))
hold on;
bode(YA(6,5));
subplot(2,2,4)
bode(YA(4,4))
hold on;
bode(YA(6,6));
title('GFL and GFM Impedance - dd axis')
legend(['GFL';'GFM'])
title('GFL and GFM Admittance - dd axis')
legend(['GFL';'GFM'])

ObjZmDss=SimplusGT.ObjSwitchInOut(ObjYmDss,14);
[~,ZAdss]= ObjZmDss.GetDSS(ObjZmDss);
figure(6);
clf;
subplot(2,2,1)
bode(ZAdss(3,3))
hold on;
bode(ZAdss(5,5));
subplot(2,2,2)
bode(ZAdss(3,4))
hold on;
bode(ZAdss(5,6));
subplot(2,2,3)
bode(ZAdss(4,3))
hold on;
bode(ZAdss(6,5));
subplot(2,2,4)
bode(ZAdss(4,4))
hold on;
bode(ZAdss(6,6));
title('GFL and GFM Impedance - dd axis')
legend(['GFL';'GFM'])

YA_tf=tf(YA);
ZA_tf=inv(YA_tf);
YAval=evalfr(YA_tf,0);
YAval=YAval(1:8,1:8);
ZAval=inv(YAval);