ObjYAss=SimplusGT.ObjDss2Ss(ObjYmDss);
[~,YA]= ObjYAss.GetSS(ObjYAss);
figure(5);
clf;
bode(YA(3,3));
grid on;
hold on;
bode(YA(5,5));
title('GFL and GFM Admittance - dd axis')
legend(['GFL';'GFM'])