% Calculate modes, residues from state-space, based on the modes and
% apparatuses the user selected.
%
% Author(s): Yue Zhu
% Modifed by: Yitong Li, Qipeng Zheng

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
    for k =1: N_Apparatus

        ZmValAll{modei}{k} = SimplusGT.Modal.ApparatusImpedanceCal(GmDSS_Cell{k}, lambda, ApparatusType{k});

        if ApparatusType{k} <= 89  % Ac apparatus
            % Using the matrix format of ResidueAll{modei}(k) should be
            % better, please change it later.
            ResidueAll{modei}{k}(1,1)=C(pout,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin);
            ResidueAll{modei}{k}(1,2)=C(pout,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin+1);
            ResidueAll{modei}{k}(2,1)=C(pout+1,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin);
            ResidueAll{modei}{k}(2,2)=C(pout+1,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin+1);           
            
        elseif ApparatusType{k} >= 1000 && ApparatusType{k} <= 1089 % Dc apparatuses
            ResidueAll{modei}{k}(1,1)=C(pout,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin);

        elseif ApparatusType{k} >= 2000 && ApparatusType{k} <= 2009 % Interlink apparatuses

            ResidueAll{modei}{k}(1,1)=C(pout,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin);
            ResidueAll{modei}{k}(1,2)=C(pout,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin+1);
            ResidueAll{modei}{k}(2,1)=C(pout+1,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin);
            ResidueAll{modei}{k}(2,2)=C(pout+1,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin+1); 

            ResidueAll{modei}{k}(1,3)=C(pout,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin+2);
            ResidueAll{modei}{k}(2,3)=C(pout+1,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin+2);
            ResidueAll{modei}{k}(3,1)=C(pout+2,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin);
            ResidueAll{modei}{k}(3,2)=C(pout+2,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin+1); 
            ResidueAll{modei}{k}(3,3)=C(pout+2,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin+2); 

        else % Floating bus and passive load: not considered           
            ResidueAll{modei}{k} = [];
        end
        pin = pin + length(ApparatusInputStr{k}); 
        pout = pout + length(ApparatusOutputStr{k});
    end
end

end