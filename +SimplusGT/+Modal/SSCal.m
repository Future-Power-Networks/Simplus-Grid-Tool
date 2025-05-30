% Calculate modes, residues from state-space, based on the modes and
% apparatuses the user selected.
%
% Author(s): Yue Zhu
% Modifed by: Yitong Li, Qipeng Zheng, George Zhang

function [Mode,ResidueAll,ZmValAll]=...
    SSCal(GminSS, N_Apparatus, ApparatusType, ModeSelect, GmDSS_Cell, ApparatusInputStr, ApparatusOutputStr)

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
    for k =1:N_Apparatus
        ZmValAll{modei}{k} = SimplusGT.Modal.ApparatusImpedanceCal(GmDSS_Cell{k}, lambda, ApparatusType{k});
        if ApparatusType{k} <= 89  % AC apparatus
            ResidueAll{modei}{k}=C(pout:pout+1,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin:pin+1);
        elseif ApparatusType{k} >= 1000 && ApparatusType{k} <= 1089 % DC apparatus
            ResidueAll{modei}{k}=C(pout,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin);
        elseif ApparatusType{k} >= 2000 && ApparatusType{k} <= 2009 % Interlink apparatus
            ResidueAll{modei}{k}=C(pout:pout+2,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin:pin+2);
        else % Floating bus and passive load: not considered           
            ResidueAll{modei}{k} = [];
        end
        pin = pin + length(ApparatusInputStr{k}); 
        pout = pout + length(ApparatusOutputStr{k});
    end
end