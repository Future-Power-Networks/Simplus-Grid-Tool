%calculate modes, residues from state-space, based on the modes and
%devicese the user selected.
%Author: Yue Zhu
function [GbMode,ResidueAll,ZmValAll,ModeTotalNum,ModeDSS,Phi_DSS, IndexSS]=...
    SSCal(GminSS, N_Bus, DeviceType, ModeSelAll, GmDSS_Cell, GsysDSS)

%% for state PF, use GsysDSS
[GsysSS, IndexSS] = SimplexPS.dss2ss(GsysDSS);
[Phi_DSS,D]=eig(GsysSS.A);
D=D/(2*pi);
ModeDSS=diag(D);

%% for impedance PF, use GminSS
A=GminSS.A;
B=GminSS.B;
C=GminSS.C;
[Phi,D]=eig(A);
D=D/(2*pi);
Psi=inv(Phi); 
GbMode=diag(D);

ModeTotalNum = length(GbMode);
ModeSelNum = length(ModeSelAll);
ResidueAll = cell(1,ModeSelNum);
ZmValAll = cell(1,ModeSelNum);

for modei=1:ModeSelNum
    ModeSel = ModeSelAll(modei);
    FreqSel = imag(GbMode(ModeSel));
    pin=1;
    pout=1;
    for k =1: N_Bus
        if DeviceType{k} <= 89  %apparatus
            %Residu calculation
            ResidueAll{modei}(k).dd=C(pout,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin);
            ResidueAll{modei}(k).dq=C(pout,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin+1);
            ResidueAll{modei}(k).qd=C(pout+1,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin);
            ResidueAll{modei}(k).qq=C(pout+1,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin+1);           
            
            %Device impedance calculation 
            GmSS_Cell = minreal(GmDSS_Cell{k});
            GmTf.dd=tf(GmSS_Cell(1,1));
            GmTf.dq=tf(GmSS_Cell(1,2));
            GmTf.qd=tf(GmSS_Cell(2,1));
            GmTf.qq=tf(GmSS_Cell(2,2));
            Gm = [GmTf.dd, GmTf.dq ; GmTf.qd, GmTf.qq];
            Zm = inv(Gm);
            %calculate the value of impedance
            1i;
            ZmValAll{modei}(k).dd = evalfr( Zm(1,1), 2*pi*FreqSel*1i);
            ZmValAll{modei}(k).dq = evalfr( Zm(1,2), 2*pi*FreqSel*1i);
            ZmValAll{modei}(k).qd = evalfr( Zm(2,1), 2*pi*FreqSel*1i);
            ZmValAll{modei}(k).qq = evalfr( Zm(2,2), 2*pi*FreqSel*1i);

            pin = pin + 4;    %4 inputs and 5 outputs.
            pout = pout + 5;
        else %floating bus: not considered.
            ResidueAll{modei}(k).dd=[];
            ZmValAll{modei}(k).dd=[];
            pin = pin + 2;    %4 inputs and 5 outputs.
            pout = pout + 2;
        end
    end
end

end