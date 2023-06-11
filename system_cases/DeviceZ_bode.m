% 
% GmDSS_Cell_base = GmDSS_Cell;
% save('GmDSS_Cell_base.mat', 'GmDSS_Cell_base');

load('GmDSS_Cell_base.mat');

k=12;
figure(3);
clf;
grid on;
P=bodeoptions;
P.Grid='on';
P.XLim={[0.1 1e4]};
P.FreqUnits='Hz';
Zmk=inv(GmDSS_Cell{k}(1:2,1:2));
Zmk_base=inv(GmDSS_Cell_base{k}(1:2,1:2));
subplot(2,2,1);
bode(Zmk(1,1),P);
hold on;
bode(Zmk_base(1,1),P);
title('impedance compare: dd axis')
legend(['Sbase=8';'Sbase=6'])

subplot(2,2,2);
bode(Zmk(1,2),P);
hold on;
bode(Zmk_base(1,2),P);

subplot(2,2,3);
bode(Zmk(2,1),P);
hold on;
bode(Zmk_base(2,1),P);

subplot(2,2,4);
bode(Zmk(2,2),P);
hold on;
bode(Zmk_base(2,2),P);

figure(4)
clf;
P=bodeoptions;
P.Grid='on';
P.XLim={[0.1 1e4]};
P.FreqUnits='Hz';
bode(GmDSS_Cell{k}(1,1),P);
hold on;
bode(GmDSS_Cell_base{k}(1,1),P);
title('admittance compare: dd axis')
legend(['Sbase=8';'Sbase=6'])