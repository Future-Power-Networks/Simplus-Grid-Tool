%% Before using. 
% State / Impedance participation factor analysis.
% Reference: Participation Analysis in Impedance Models: The Grey-Box Approach for Power System Stability
%  
% Author(s): Yue Zhu

%% Notes
% Before running this program, you need to config the analysis from
% ModalConfig.xlsx, which should be located in the toolbox root folder.
% Note1: device numbering keeps the same as bus numbering. For example: the device
% on bus7 will always be named as Device7.
% Note2: The final results will be saved in MdLayer1, MdLayer2, MdLayer3,
% MdMode, MdStatePF.

%% Basic
%Basic infomation acquirement.
FileModal=['ModalConfig_' Name_Netlist '.xlsx'];
clear MdLayer1;
clear MdLayer2;
clear MdLayer3;
clear MdStatePF;
clear MdMode;
%tic
%read Modal config file.
[AxisSel, DeviceSelL12, ModeSelAll, DeviceSelL3All,StateSel_DSS, ModeSel_DSS] = ...
    SimplexPS.Modal.ExcelRead(FileModal, N_Bus, DeviceType, GminSS);
[StatePFEnable, BodeEnable, Layer12Enable, Layer3Enable] = ...
    SimplexPS.Modal.EnablingRead(FileModal); %Enablling control.

%check for illegal selection.
SimplexPS.Modal.DataCheck(AxisSel, DeviceSelL12, ModeSelAll, DeviceSelL3All,...
    StateSel_DSS, ModeSel_DSS,BodeEnable,Layer12Enable,Layer3Enable,StatePFEnable);

ModeSelNum = length(ModeSelAll);
%get ResidueAll, ZmValAll.
[MdMode,ResidueAll,ZmValAll,ModeTotalNum,ModeDSS,Phi_DSS, IndexSS]=...
    SimplexPS.Modal.SSCal(GminSS, N_Bus, DeviceType, ModeSelAll, GmDSS_Cell, GsysDSS, DeviceInputStr, DeviceOutputStr);

%% Impedance Participation Factor
%Analysis.
if BodeEnable ==1
    fprintf('plotting bode diagram for selected whole-system admittance...\n')
    SimplexPS.Modal.BodeDraw(DeviceSelL12, AxisSel, GminSS, DeviceType, N_Bus, DeviceInputStr, DeviceOutputStr);
end

for modei=1:ModeSelNum
    Residue = ResidueAll{modei};
    ZmVal = ZmValAll{modei};
    FreqSel = imag(MdMode(ModeSelAll(modei)));
    if Layer12Enable ==1
        fprintf('Calculating Modal Analysis Layer1&2 and plotting the results...\n')
        [Layer1, Layer2] = SimplexPS.Modal.MdLayer12(Residue,ZmVal,N_Bus,...
            DeviceType,modei,DeviceSelL12,FreqSel,MdMode(ModeSelAll(modei)));
        MdLayer1(modei).mode = [num2str(FreqSel),'~Hz'];
        MdLayer2(modei).mode = [num2str(FreqSel),'~Hz'];
        for count = 1: length(DeviceSelL12)
            MdLayer1(modei).result(count).Device={['Device',num2str(DeviceSelL12(count))]};
            MdLayer1(modei).result(count).Abs_Max=Layer1(count);
            MdLayer2(modei).result(count).Device={['Device',num2str(DeviceSelL12(count))]};
            MdLayer2(modei).result(count).DeltaLambdaReal=Layer2.real(count);
            MdLayer2(modei).result(count).DeltaLambdaImag=Layer2.imag(count);
            MdLayer2(modei).result(count).DeltaLambdaRealpu=Layer2.real_pu(count);
            MdLayer2(modei).result(count).DeltaLambdaImagpu=Layer2.imag_pu(count);
        end
    end
    if Layer3Enable ==1
        fprintf('Calculating Modal Layer3...\n')
        MdLayer3(modei).mode = [num2str(FreqSel),'~Hz'];
        MdLayer3(modei).result = SimplexPS.Modal.MdLayer3(Residue,ZmVal,...
        FreqSel,DeviceType,DeviceSelL3All,Para,PowerFlow,Ts);
    end
end

%% State Participation Factor
if StatePFEnable == 1
%Phi_DSS;
%StateSel_DSS;
DeviceStateTotal = 0;
for Di=1:length(DeviceStateStr)
    DeviceStateTotal = DeviceStateTotal + length(DeviceStateStr{Di});
end

for modei = 1: length(ModeSel_DSS)
    FreqSel = imag(ModeDSS(ModeSel_DSS(modei)));
    ModeSel = ModeSel_DSS(modei);
    MdStatePF(modei).mode = [num2str(FreqSel),'~Hz'];
    for statei=1:length(StateSel_DSS)
            StateSel = IndexSS(StateSel_DSS(statei)); % this is used to print the correct name
            if StateSel > DeviceStateTotal %belong to line.
                StateName = {'Line'};
            else %belong to device
                StateSel_ = StateSel;
                for Di = 1:N_Device
                    if StateSel_ <= length(DeviceStateStr{Di})
                        StateName = {['Device',num2str(Di)]};
                        break;
                    else
                        StateSel_ = StateSel_ - length(DeviceStateStr{Di});
                    end
                end
            end
            MdStatePF(modei).result(statei).Device = StateName;
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
toc
