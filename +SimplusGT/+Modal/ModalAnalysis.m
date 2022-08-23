%% Before using. 
% State / Impedance participation factor analysis.
% Reference: Participation Analysis in Impedance Models: The Grey-Box Approach for Power System Stability
%  
% Author(s): Yue Zhu

%% Notes
% Before running this program, you need to config the analysis from
% ModalConfig.xlsx, which should be located in the toolbox root folder.
% Note1: apparatus numbering keeps the same as bus numbering. For example: the apparatus
% on bus7 will always be named as Apparatus7.
% Note2: The final results will be saved in MdLayer1, MdLayer2, MdLayer3,
% MdMode, MdStatePF.

%% Basic
% Change suffix
UserData_Modal = UserDataName;
UserData_Modal = strrep(UserData_Modal,'.xlsm','');
UserData_Modal = strrep(UserData_Modal,'.xlsx','');
UserData_Modal = strrep(UserData_Modal,'.xls','');
UserData_Modal = strrep(UserData_Modal,'.json','');
UserData_Modal = [UserData_Modal,'.xlsx'];

%Basic infomation acquirement.
FileModal=[cd '\Examples\ParticipationAnalysis\ModalConfig_' UserData_Modal];
clear MdLayer1;
clear MdLayer2;
clear MdLayer3;
clear MdStatePF;
clear MdMode;

%read Modal config file.
[AxisSel, ApparatusSelL12, ModeSelAll, ApparatusSelL3All,StateSel_DSS, ModeSel_DSS] = ...
    SimplusGT.Modal.ExcelRead(FileModal, NumBus, ApparatusType, GsysSs);
[StatePFEnable, BodeEnable, Layer12Enable, Layer3Enable] = ...
    SimplusGT.Modal.EnablingRead(FileModal); %Enablling control.

%check for illegal selection.
SimplusGT.Modal.DataCheck(AxisSel, ApparatusSelL12, ModeSelAll, ApparatusSelL3All,...
    StateSel_DSS, ModeSel_DSS,BodeEnable,Layer12Enable,Layer3Enable,StatePFEnable);

ModeSelNum = length(ModeSelAll);
%get ResidueAll, ZmValAll.
[MdMode,ResidueAll,ZmValAll,ModeTotalNum,ModeDSS,Phi_DSS, IndexSS]=...
    SimplusGT.Modal.SSCal(GsysSs, NumBus, ApparatusType, ModeSelAll, GmDssCell, GsysDss, ApparatusInputStr, ApparatusOutputStr);

%% Impedance Participation Factor
%Analysis.
if BodeEnable ==1
    fprintf('plotting bode diagram for selected whole-system admittance...\n')
    SimplusGT.Modal.BodeDraw(ApparatusSelL12, AxisSel, GsysSs, ApparatusType, ApparatusBus, NumBus, ApparatusInputStr, ApparatusOutputStr);
end

for modei=1:ModeSelNum
    Residue = ResidueAll{modei};
    ZmVal = ZmValAll{modei};
    FreqSel = imag(MdMode(ModeSelAll(modei)));
    if Layer12Enable ==1
        fprintf('Calculating Modal Analysis Layer1&2 and plotting the results...\n')
        [Layer1, Layer2] = SimplusGT.Modal.MdLayer12(Residue,ZmVal,N_Bus,...
            ApparatusType,modei,ApparatusSelL12,FreqSel,MdMode(ModeSelAll(modei)));
        MdLayer1(modei).mode = [num2str(FreqSel),'~Hz'];
        MdLayer2(modei).mode = [num2str(FreqSel),'~Hz'];
        for count = 1: length(ApparatusSelL12)
            MdLayer1(modei).result(count).Apparatus={['Apparatus',num2str(ApparatusSelL12(count))]};
            MdLayer1(modei).result(count).Abs_Max=Layer1(count);
            MdLayer2(modei).result(count).Apparatus={['Apparatus',num2str(ApparatusSelL12(count))]};
            MdLayer2(modei).result(count).DeltaLambdaReal=Layer2.real(count);
            MdLayer2(modei).result(count).DeltaLambdaImag=Layer2.imag(count);
            MdLayer2(modei).result(count).DeltaLambdaRealpu=Layer2.real_pu(count);
            MdLayer2(modei).result(count).DeltaLambdaImagpu=Layer2.imag_pu(count);
        end
    end
    if Layer3Enable ==1
        fprintf('Calculating Modal Layer3...\n')
        MdLayer3(modei).mode = [num2str(FreqSel),'~Hz'];
        MdLayer3(modei).result = SimplusGT.Modal.MdLayer3(Residue,ZmVal,FreqSel,ApparatusType,...
                ApparatusSelL3All,Para,PowerFlow,Ts,ApparatusBus,ListBus);
    end
end

%% State Participation Factor
if StatePFEnable == 1
%Phi_DSS;
%StateSel_DSS;
ApparatusStateTotal = 0;
for Di=1:length(ApparatusStateStr)
    ApparatusStateTotal = ApparatusStateTotal + length(ApparatusStateStr{Di});
end

for modei = 1: length(ModeSel_DSS)
    FreqSel = imag(ModeDSS(ModeSel_DSS(modei)));
    ModeSel = ModeSel_DSS(modei);
    MdStatePF(modei).mode = [num2str(FreqSel),'~Hz'];
    for statei=1:length(StateSel_DSS)
            StateSel = IndexSS(StateSel_DSS(statei)); % this is used to print the correct name
            if StateSel > ApparatusStateTotal %belong to line.
                StateName = {'Line'};
            else %belong to apparatus
                StateSel_ = StateSel;
                for Di = 1:NumApparatus
                    if StateSel_ <= length(ApparatusStateStr{Di})
                        StateName = {['Apparatus',num2str(Di)]};
                        break;
                    else
                        StateSel_ = StateSel_ - length(ApparatusStateStr{Di});
                    end
                end
            end
            MdStatePF(modei).result(statei).Apparatus = StateName;
            MdStatePF(modei).result(statei).State = SysStateString(StateSel);
            % after printing the correct name, change StateSel back to match with GsysSS.
            StateSel = StateSel_DSS(statei); 
            Psi_DSS = inv(Phi_DSS);
            StatePF = Phi_DSS(StateSel,ModeSel) * Psi_DSS(ModeSel,StateSel);            
            MdStatePF(modei).result(statei).PF = StatePF;
            MdStatePF(modei).result(statei).PF_ABS = abs(StatePF);
            MdStatePF(modei).result(statei).PF_Real = real(StatePF);
            MdStatePF(modei).result(statei).PF_Imag = imag(StatePF);
    end
    
end
end

