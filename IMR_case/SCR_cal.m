Z14=0.01+0.1*1i; Y14=Z14^(-1);
Z24=0.02+0.2*1i; Y24=Z24^(-1);
Z34=0.02+0.2*1i; Y34=Z34^(-1);
Zc2=1/(0.3*1i); Yc2=Zc2^(-1);
Zc3=1/(0.3*1i); Yc3=Zc3^(-1);


Z1 = (Z34+Zc3)*Z14/(Z34+Zc3+Z14);
Z2 = Z1+Z24;
Zth2=Z2*Zc2/(Z2+Zc2);

Z3 = (Z24+Zc2)*Z14/(Z24+Zc2+Z14);
Z4 = Z3+Z34;
Zth3=Z4*Zc2/(Z4+Zc2);

SCR2 = 1/abs(Zth2)/0.32
SCR2 = 1/abs(Zth3)/0.32


Ynodal = [Y14,  0,          0,  -Y14;...
          0,    Yc2+Y24,    0,  -Y24;...
          0,    0,          Yc3+Y34,    -Y34;...
          -Y14, -Y24,       -Y34,   1e-6];
Ynodal_red = Ynodal(1:3,1:3)-Ynodal(1:3,4)*Ynodal(4,4)^(-1)*Ynodal(4,1:3)

Znodal = inv(Ynodal_red);

IF32= abs(Znodal(3,2)/Znodal(2,2));
IF23 = abs(Znodal(2,3)/Znodal(3,3));
ESCR2 = (1/abs(Zth2))/(0.32+0.32*IF32)
ESCR3 = (1/abs(Zth3))/(0.32+0.32*IF23)
