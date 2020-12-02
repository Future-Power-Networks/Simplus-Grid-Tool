function ZmVal = DeviceImpedanceCal(GmDSS_Cell, FreqSel, DeviceSel)
GmSS = minreal(GmDSS_Cell);
GmTf.dd=tf(GmSS(1,1));
GmTf.dq=tf(GmSS(1,2));
GmTf.qd=tf(GmSS(2,1));
GmTf.qq=tf(GmSS(2,2));
Gm = [GmTf.dd, GmTf.dq ; GmTf.qd, GmTf.qq];
Zm = inv(Gm);
%calculate the value of impedance
1i;
ZmVal.dd = evalfr( Zm(1,1), 2*pi*FreqSel*1i);
ZmVal.dq = evalfr( Zm(1,2), 2*pi*FreqSel*1i);
ZmVal.qd = evalfr( Zm(2,1), 2*pi*FreqSel*1i);
ZmVal.qq = evalfr( Zm(2,2), 2*pi*FreqSel*1i);
end