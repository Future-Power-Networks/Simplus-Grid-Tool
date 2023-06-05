function [MdLayer1, MdLayer2, MdLayer3, MdStatePF, MdMode, MdSensResult, ZminSS]=ModalAnalysisExe(UserdataModal)

%% invoke variables from original workspace
N_Apparatus = evalin('base', 'N_Apparatus');
N_Bus = evalin('base', 'N_Bus');
ApparatusType = evalin('base', 'ApparatusType');
ApparatusBus = evalin('base', 'ApparatusBus');
ApparatusInputStr = evalin('base', 'ApparatusInputStr');
ApparatusOutputStr = evalin('base', 'ApparatusOutputStr');
ApparatusStateStr=evalin('base', 'ApparatusStateStr');
SysStateString=evalin('base', 'SysStateString');
GminSS = evalin('base', 'GminSS');
GmDSS_Cell = evalin('base', 'GmDSS_Cell');
GsysDSS = evalin('base', 'GsysDSS');
Para = evalin('base', 'Para');
ApparatusPowerFlow = evalin('base', 'ApparatusPowerFlow');
Ts = evalin('base', 'Ts');
ListBus = evalin('base', 'ListBus');

GmObj=evalin('base', 'GmObj');
YbusObj=evalin('base', 'YbusObj');
Port_i=evalin('base', 'Port_i');
Port_v=evalin('base', 'Port_v');


CaseStudy=evalin('base', 'CaseStudy'); % for case study selection

%% read Modal config file.

[AxisSel, ApparatusSelL12, ModeSelAll, ApparatusSelL3All,StateSel_DSS, ModeSel_DSS] = ...
    SimplusGT.Modal.ExcelRead(UserdataModal, N_Apparatus, ApparatusType, GminSS);

[StatePFEnable, BodeEnable, Layer12Enable, Layer3Enable, SensEnable] = ...
    SimplusGT.Modal.EnablingRead(UserdataModal); %Enablling control.

%check for illegal selection.
SimplusGT.Modal.DataCheck(AxisSel, ApparatusSelL12, ModeSelAll, ApparatusSelL3All,...
    StateSel_DSS, ModeSel_DSS,BodeEnable,Layer12Enable,Layer3Enable,StatePFEnable);

ModeSelNum = length(ModeSelAll);
%get ResidueAll, ZmValAll.
[MdMode,ResidueAll,ZmValAll,~,ModeDSS,Phi_DSS, IndexSS]=...
    SimplusGT.Modal.SSCal(GminSS, N_Apparatus, ApparatusType, ModeSelAll, GmDSS_Cell, GsysDSS, ApparatusInputStr, ApparatusOutputStr);

%% Impedance Participation Factor
%Analysis.
if BodeEnable ==1
    fprintf('plotting bode diagram for selected whole-system admittance...\n')
    SimplusGT.Modal.BodeDraw(ApparatusSelL12, AxisSel, GminSS, ApparatusType, ApparatusBus, N_Apparatus, ApparatusInputStr, ApparatusOutputStr);
end

for modei=1:ModeSelNum
    Residue = ResidueAll{modei};
    ZmVal = ZmValAll{modei};
    FreqSel = imag(MdMode(ModeSelAll(modei)));
    if Layer12Enable ==1
        fprintf('Calculating Modal Analysis Layer1&2 and plotting the results...\n')
        [Layer1, Layer2] = SimplusGT.Modal.MdLayer12(Residue,ZmVal,N_Apparatus,ApparatusBus,...
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
    else % Layer-1-2 not selected
        MdLayer2=0;
    end
    if Layer3Enable ==1
        fprintf('Calculating Modal Layer3...\n')
        MdLayer3(modei).mode = [num2str(FreqSel),'~Hz'];
        Mode_Hz = MdMode(ModeSelAll(modei));
        MdLayer3(modei).result = SimplusGT.Modal.MdLayer3(Residue,ZmVal,Mode_Hz,ApparatusType,...
                ApparatusSelL3All,Para,ApparatusPowerFlow,Ts,ApparatusBus,ListBus);
    else % Layer3 not Enabled.
        MdLayer3 = 0;
    end
end

%% State Participation Factor
if StatePFEnable == 1
    fprintf('Calculating state-space participation factor...\n')
%Phi_DSS;
%StateSel_DSS;
ApparatusStateTotal = 0;
for Di=1:length(ApparatusStateStr)
    ApparatusStateTotal = ApparatusStateTotal + length(ApparatusStateStr{Di});
end

Psi_DSS = inv(Phi_DSS);
for modei = 1: length(ModeSel_DSS)
    FreqSel = imag(ModeDSS(ModeSel_DSS(modei)));
    ModeSel = ModeSel_DSS(modei);
    MdStatePF(modei).mode = [num2str(ModeDSS(ModeSel_DSS(modei))),'~Hz'];
    for statei=1:length(StateSel_DSS)
            StateSel = IndexSS(StateSel_DSS(statei)); % this is used to print the correct name
            if StateSel > ApparatusStateTotal %belong to line.
                StateName = {'Line'};
            else %belong to apparatus
                StateSel_ = StateSel;
                for Di = 1:N_Apparatus
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
            StatePF = Phi_DSS(StateSel,ModeSel) * Psi_DSS(ModeSel,StateSel);            
            MdStatePF(modei).result(statei).PF = StatePF;
            MdStatePF(modei).result(statei).PF_ABS = abs(StatePF);
            MdStatePF(modei).result(statei).PF_Real = real(StatePF);
            MdStatePF(modei).result(statei).PF_Imag = imag(StatePF);
    end
    
end
else % State participation not enabled.
    MdStatePF = 0;
end


%% Admittance Sensitivity Analysis (node + branch) = -Residues matrix
% Ybus: nodal addmittance matrix without apparatus
% Ynodal: nodal admittance matrix with apparatus
% Yre: a rearranged admittance: diagonal---node element admittance;
%                               off-diagonal------branch element admittance
MdSensResult=[];

if SensEnable == 1 % if enable
fprintf('Calculating whole-system eigenvalue sensitivity...\n')
ZminSS = SimplusGT.WholeSysZ_cal(GmObj,YbusObj,Port_i, Port_v);
[~,D]=eig(ZminSS.A);
ZMode_rad=diag(D);
ZMode_Hz=ZMode_rad/2/pi;
for modei=1:ModeSelNum
    
    % modes of Zsys may not completely matched with Ysys (a small error due to float calculation), so we need to
    % find the same (neareast) mode in Zsys first.
    m_tol = 1e-3; % tolerance
     
    lambda_sel = MdMode(ModeSelAll(modei));
    Ek = 0;
    for n = 1: length(ZMode_Hz) % Mode calculated from Y may in different order from calculated from Z
        if abs(lambda_sel-ZMode_Hz(n))<=m_tol
            Ek = n;
        end
        if n == length(ZMode_Hz) && abs(lambda_sel-ZMode_Hz(n))>m_tol && Ek == 0 % not matched
            error('eigenvalue mismatch');
        end
    end
    if Ek==0
        error('eigenvalue mismatch');
    end
    FreqSel = imag(ZMode_Hz(Ek));
    Mode_rad = ZMode_rad(Ek);
   
    fprintf('Calculating sensitivity matrix...\n')
    [SensMatrix, ~, Ynodal_val, Yre_val, SensMat_exp] ...
        =SimplusGT.Modal.SensitivityCal(ZminSS,Ek,Mode_rad,YbusObj); % get the sensitivity matrix and some admittance matrices.        
    
    fprintf('Calculating sensitivity Layer-2...\n')
    [SensLayer1_val, SensLayer2_val,Layer12] = SimplusGT.Modal.SensLayer12(SensMatrix,Yre_val,modei,ZMode_Hz(Ek));
    fprintf('Calculating sensitivity Layer-3...\n')
    Line_sel = [15];% See ListLineNew!.
    if CaseStudy>=4 && CaseStudy<=6
        Line_sel = [37];
    elseif CaseStudy>6
        Line_sel = [2];
    end
    [SensLayer3_app,SensLayer3_bus] = SimplusGT.Modal.SensLayer3(SensMatrix,Mode_rad,ApparatusSelL3All,Line_sel);
    
    % small-signal strength
    ZmVal = ZmValAll{modei};
    ResidueY = ResidueAll{modei};
    AMR=zeros(N_Bus,1);
    IMR=zeros(N_Bus,1);
    IMRx=zeros(N_Bus,1);
    Z_layer2=zeros(N_Bus,1);
    SigmaMag= abs(real(Mode_rad));
    for k=1:N_Bus
        if ApparatusType{k} ~= 100 % not a floating bus
            Res_k = -SensMat_exp(2*k-1:2*k, 2*k-1:2*k);
            YA_k=evalfr(GmDSS_Cell{k}(1:2,1:2),Mode_rad);
            AMR(k)=SigmaMag / abs(-1*( Res_k(1,1)*YA_k(1,1) + Res_k(1,2)*YA_k(2,1) + Res_k(2,1)*YA_k(1,2) + Res_k(2,2)*YA_k(2,2) ));
            %AMR(k)= 1+SigmaMag / (norm(Res_k,'fro')*norm(YA_k,'fro'));
            %IMR(k) = SigmaMag / (SimplusGT.Frobenius_norm_dq(ResidueY(k)) * SimplusGT.Frobenius_norm_dq(ZmVal(k)));
            %IMRx(k) = SigmaMag / (abs(-1*SimplusGT.inner_product_dq(ResidueY(k),ZmVal(k))));
            IMR(k)= SigmaMag/abs(-1*SimplusGT.inner_product_dq(ResidueY(k),ZmVal(k)));    
            %Y_layer2(k)= -1*( Res_k(1,1)*YA_k(1,1) + Res_k(1,2)*YA_k(2,1) + Res_k(2,1)*YA_k(1,2) + Res_k(2,2)*YA_k(2,2) );
            %Z_layer2(k)= -1*SimplusGT.inner_product_dq(ResidueY(k),ZmVal(k));
            %SimplusGT.Frobenius_norm_dq(SensMatrix(1,1))
        end
    end
    MdSensResult(modei).mode = [num2str(FreqSel),'~Hz'];
    MdSensResult(modei).SensMatrix = SensMatrix;
    MdSensResult(modei).Ynodal_val = Ynodal_val;
    MdSensResult(modei).Yre_val = Yre_val;
    MdSensResult(modei).Layer1_mat = SensLayer1_val;
    MdSensResult(modei).Layer2_mat = SensLayer2_val;
    MdSensResult(modei).Layer12 = Layer12;
    MdSensResult(modei).Layer3_app = SensLayer3_app;
    MdSensResult(modei).Layer3_bus = SensLayer3_bus;
    MdSensResult(modei).AMR = AMR;
    MdSensResult(modei).IMR= IMR;
    MdSensResult(modei).IMRx=IMRx;
    %MdSensResult(modei).Z_Layer2=Z_layer2;
    %MdSensResult.Sen
%     figure(modei+100);
%     SensLayer1_val(find(SensLayer1_val==0))=nan;
%     heatmap(SensLayer1_val);
%     title('Heatmap of Grey-box Layer-1')
%     figure(modei+101);
%     SensLayer2_val(find(SensLayer2_val==0))=nan;
%     heatmap(real(SensLayer2_val));
%     title('Heatmap of Grey-box Layer-2 real part (damping)')
    
end
end
% Sensitivity Layer 1


end % end of function