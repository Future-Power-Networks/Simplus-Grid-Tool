% 1) in 4bus case, pole map of GsysNew and GminSS are matched, while det()
% function turns out to be wrong.
% 2) in 14bus case, pole map of GsysNew and GminSS are not matched. 
% To change the case into 14bus: in SimplexPS/Toolbox/Main.m, change 'CustomerData.xlsx' into
% 'CustomerData2.xlsx', which is the 14 bus case.

clear GsysSS;
clear IndexSS;
[GsysSS, IndexSS] = SimplexPS.dss2ss(GsysDSS);

%D=eig(GsysNew.A)/(2*pi);
D=eig(GsysSS)/(2*pi);
%GsysNew2=minreal(GsysNew);
%E = eig(GsysNew2.A)/(2*pi);

figure(1001)
hold on
subplot(1,2,1)
scatter(real(D),imag(D),'x','LineWidth',1.5); hold on; grid on;
