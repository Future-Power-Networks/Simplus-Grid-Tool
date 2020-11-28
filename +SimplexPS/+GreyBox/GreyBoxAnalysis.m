%Greybox approach for participation analysis.
% Author(s): Yue Zhu

%% Notes
% Before running this program, you need to config the analysis from
% GreyboxConfig.xlsx, which should be located in the toolbox root folder.
% Note1: device numbering keeps the same as bus numbering. For example: the device
% on bus7 will always be named as Device7.
% Note2: The final results will be saved in GbLayer1, GbLayer2, GbLayer3,
% GbMode and GbResidue.

%% Basic
%Basic infomation acquirement.
filename='GreyBoxConfig.xlsx';
BodeEnable = 1;
Layer12Enable = 1;
Layer3Enable = 1;

[AxisSel, DeviceSel12, ModeSel, DeviceSel3,StateSel_DSS, ModeSel_DSS] = ...
    SimplexPS.GreyBox.ExcelRead(filename, N_Bus, DeviceType, GminSS);
ModeSelNum = length(ModeSelAll);
%get ResidueAll, ZmValAll.
[GbMode,ResidueAll,ZmValAll,ModeTotalNum,ModeDSS,Phi_DSS]=...
    SimplexPS.GreyBox.SSCal(GminSS, N_Bus, DeviceType, ModeSelAll, GmDSS_Cell, GsysDSS);

%% Impedance Participation Factor
%Analysis.
if BodeEnable ==1
    fprintf('plotting bode diagram for selected whole-system admittance...\n')
    SimplexPS.GreyBox.BodeDraw(DeviceSelL12, AxisSel, GminSS, DeviceType, N_Bus);
end

clear GbLayer1;
clear GbLayer2;
clear GbLayer3;
for modei=1:ModeSelNum
    Residue = ResidueAll{modei};
    ZmVal = ZmValAll{modei};
    FreqSel = imag(GbMode(ModeSelAll(modei)));
    if Layer12Enable ==1
        fprintf('Calculating GreyBox Layer1&2 and plotting the results...\n')
        [Layer1, Layer2] = SimplexPS.GreyBox.GbLayer12(Residue,ZmVal,N_Bus,...
            DeviceType,modei,DeviceSelL12,FreqSel);
        GbLayer1(modei).mode = [num2str(FreqSel),'~Hz'];
        GbLayer2(modei).mode = [num2str(FreqSel),'~Hz'];
        for count = 1: length(DeviceSelL12)
            GbLayer1(modei).result(count).Device={['Device',num2str(DeviceSelL12(count))]};
            GbLayer1(modei).result(count).Abs_Max=Layer1(count);
            GbLayer2(modei).result(count).Device={['Device',num2str(DeviceSelL12(count))]};
            GbLayer2(modei).result(count).DeltaLambdaReal=Layer2.real(count);
            GbLayer2(modei).result(count).DeltaLambdaImag=Layer2.imag(count);
        end
    end
    if Layer3Enable ==1
        fprintf('Calculating GreyBox Layer3...\n')
        GbLayer3(modei).mode = [num2str(FreqSel),'~Hz'];
        GbLayer3(modei).result = SimplexPS.GreyBox.GbLayer3(Residue,ZmVal,...
        FreqSel,DeviceType,DeviceSelL3All,Para,PowerFlow,Ts);
    end
end

%% State Participation Factor
Phi_DSS;
StateSel_DSS;
DeviceStateTotal = 0;
for Di=1:length(DeviceStateStr)
    DeviceStateTotal = DeviceStateTotal + length(DeviceStateStr{Di});
end

for modei = 1: length(ModeSel_DSS)
    FreqSel = imag(ModeDSS(ModeSel_DSS(modei)));
    ModeSel = ModeSel_DSS(modei);
    GbStatePF(modei).mode = [num2str(FreqSel),'~Hz'];
    for statei=1:length(StateSel_DSS)
            StateSel = StateSel_DSS(statei);
            DeviceNum = length(DeviceStateStr);
            if StateSel > DeviceStateTotal %belong to line.
                StateName = {'Line'};
            else %belong to device
                StateSel_ = StateSel;
                for Di = 1:DeviceNum
                    if StateSel_ <= length(DeviceStateStr{Di})
                        StateName = {['Device',num2str(Di)]};
                        break;
                    else
                        StateSel_ = StateSel_ - length(DeviceStateStr{Di});
                    end
                end
            end
            StatePF = Phi_DSS(ModeSel,StateSel) * Phi_DSS(StateSel,ModeSel);
            GbStatePF(modei).result(statei).Device = StateName;
            GbStatePF(modei).result(statei).State = SysStateString(StateSel);
            GbStatePF(modei).result(statei).PF = StatePF;
            GbStatePF(modei).result(statei).PF_ABS = abs(StatePF);
    end
    

end
  

