%calculate modes, residues from state-space, based on the modes and
%apparatuses the user selected.
%Author: Yue Zhu
function [Mode,ResidueAll,ZmValAll]=...
    SSCal(GminSS, N_Apparatus, ApparatusType, ModeSelect, GmDSS_Cell, ApparatusInputStr, ApparatusOutputStr)

%% for state PF, use GsysDSS
%[GsysSS, IndexSS] = SimplusGT.dss2ss(GsysDSS);
%[Phi_DSS,D]=eig(GsysSS.A);
%D=D/(2*pi);
%ModeDSS=diag(D);

%% for impedance PF, use GminSS
A=GminSS.A;
B=GminSS.B;
C=GminSS.C;
[Phi,D]=eig(A);
D=D/(2*pi);
Psi=inv(Phi); 
Mode=diag(D);

ModeSelNum = length(ModeSelect);
ResidueAll = cell(1,ModeSelNum);
ZmValAll = cell(1,ModeSelNum);

for modei=1:ModeSelNum
    ModeSel = ModeSelect(modei);
    lambda = Mode(ModeSel)*2*pi;
    pin=1;
    pout=1;
    for k =1: N_Apparatus
        if ApparatusType{k} <= 89  %apparatus
            %Residu calculation
            ResidueAll{modei}(k).dd=C(pout,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin);
            ResidueAll{modei}(k).dq=C(pout,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin+1);
            ResidueAll{modei}(k).qd=C(pout+1,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin);
            ResidueAll{modei}(k).qq=C(pout+1,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin+1);           
            
            1i;
            GmTf.dd=evalfr(GmDSS_Cell{k}(1,1),lambda);
            GmTf.dq=evalfr(GmDSS_Cell{k}(1,2),lambda);
            GmTf.qd=evalfr(GmDSS_Cell{k}(2,1),lambda);
            GmTf.qq=evalfr(GmDSS_Cell{k}(2,2),lambda);
            Gm = [GmTf.dd, GmTf.dq ; GmTf.qd, GmTf.qq];
            Zm = inv(Gm);
            ZmValAll{modei}(k).dd = Zm(1,1);
            ZmValAll{modei}(k).dq = Zm(1,2);
            ZmValAll{modei}(k).qd = Zm(2,1);
            ZmValAll{modei}(k).qq = Zm(2,2);
        else %floating bus and passive load: not considered           
            ResidueAll{modei}(k).dd=[];
            ZmValAll{modei}(k).dd=[];
        end
        pin = pin + length(ApparatusInputStr{k});    %4 inputs and 5 outputs.
        pout = pout + length(ApparatusOutputStr{k});
    end
end

end