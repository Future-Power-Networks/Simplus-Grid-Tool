function ZmVal = ApparatusImpedanceCal(GmDSS_Cell, FreqSel, ApparatusSel)

1i;
GmTf.dd=evalfr(GmDSS_Cell(1,1),2*pi*FreqSel*1i);
GmTf.dq=evalfr(GmDSS_Cell(1,2),2*pi*FreqSel*1i);
GmTf.qd=evalfr(GmDSS_Cell(2,1),2*pi*FreqSel*1i);
GmTf.qq=evalfr(GmDSS_Cell(2,2),2*pi*FreqSel*1i);

Gm = [GmTf.dd, GmTf.dq ; GmTf.qd, GmTf.qq];
Zm = inv(Gm);
%calculate the value of impedance
ZmVal.dd = Zm(1,1);
ZmVal.dq = Zm(1,2);
ZmVal.qd = Zm(2,1);
ZmVal.qq = Zm(2,2);
end