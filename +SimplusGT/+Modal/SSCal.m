%calculate modes, residues from state-space, based on the modes and
%apparatuses the user selected.
%Author: Yue Zhu
function [GbMode,ResidueAll,ZmValAll,ModeTotalNum,ModeDSS,Phi_DSS, IndexSS]=...
    SSCal(GminSS, N_Bus, ApparatusType, ModeSelAll, GmDSS_Cell, GsysDSS, ApparatusInputStr, ApparatusOutputStr)

%% for state PF, use GsysDSS
[GsysSS, IndexSS] = SimplusGT.dss2ss(GsysDSS);
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
        if ApparatusType{k} <= 89  %apparatus
            %Residu calculation
            ResidueAll{modei}(k).dd=C(pout,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin);
            ResidueAll{modei}(k).dq=C(pout,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin+1);
            ResidueAll{modei}(k).qd=C(pout+1,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin);
            ResidueAll{modei}(k).qq=C(pout+1,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin+1);           
            
            1i;
            GmTf.dd=evalfr(GmDSS_Cell{k}(1,1),2*pi*FreqSel*1i);
            GmTf.dq=evalfr(GmDSS_Cell{k}(1,2),2*pi*FreqSel*1i);
            GmTf.qd=evalfr(GmDSS_Cell{k}(2,1),2*pi*FreqSel*1i);
            GmTf.qq=evalfr(GmDSS_Cell{k}(2,2),2*pi*FreqSel*1i);
            Gm = [GmTf.dd, GmTf.dq ; GmTf.qd, GmTf.qq];
            Zm = inv(Gm);
            ZmValAll{modei}(k).dd = Zm(1,1);
            ZmValAll{modei}(k).dq = Zm(1,2);
            ZmValAll{modei}(k).qd = Zm(2,1);
            ZmValAll{modei}(k).qq = Zm(2,2);
            pin = pin + length(ApparatusInputStr{k});    %4 inputs and 5 outputs.
            pout = pout + length(ApparatusOutputStr{k});
            
        else %floating bus and passive load: not considered           
            ResidueAll{modei}(k).dd=[];
            ZmValAll{modei}(k).dd=[];
            pin = pin + length(ApparatusInputStr{k});
            pout = pout + length(ApparatusOutputStr{k});
        end
    end
end

end