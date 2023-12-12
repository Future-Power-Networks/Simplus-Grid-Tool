
N_Apparatus = evalin('base', 'N_Apparatus');
N_Bus = evalin('base', 'N_Bus');
ApparatusType = evalin('base', 'ApparatusType');
ApparatusBus = evalin('base', 'ApparatusBus');
ApparatusInputStr = evalin('base', 'ApparatusInputStr');
ApparatusOutputStr = evalin('base', 'ApparatusOutputStr');
ApparatusStateStr=evalin('base', 'ApparatusStateStr');
GminStateStr=evalin('base', 'GminStateStr');
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


A=GminSS.A;
[~,D]=eig(A);
D_Hz=diag(D/(2*pi));

m=1;
ModeSelect=0;
for i=1:length(D_Hz)
    if imag(D_Hz(i))<0 
        %negative imaginary part - ignored
    elseif imag(D_Hz(i))<1
        % modes smaller than 1Hz are not included
    elseif abs(imag(D_Hz(i))-Fbase)<1
        % modes close to base frequency are not included
    %elseif abs(real(D_Hz(i)))>20
        % with a real part that is very large
    elseif imag(D_Hz(i))>500
        % larger than 500 Hz are out of consideration
    %elseif abs(imag(D_Hz(i))) < (abs(real(D_Hz(i))) * 7)
        % larger than 30% damping ratio
    else
        ModeSelect(m)=i;
        m=m+1;
    end
end

[MdMode,ResidueAll,ZmValAll,~,~,~, ~]=...
    SimplusGT.Modal.SSCal(GminSS, N_Apparatus, ApparatusType, ModeSelect, GmDSS_Cell, GsysDSS, ApparatusInputStr, ApparatusOutputStr);

%% config C-IMR: Floating bus with maximum IMR value, SG with no effect on heat map, and IBR has the critical effect.
j=1;
clear CIMR CIMR2;
for k=1:N_Bus
    if ApparatusType{k}<30 || ApparatusType{k}>=40
        CIMR2(j).device = k;
        CIMR2(j).value = log10(100);
        CIMR2(j).mode = 0;
        j=j+1;
    end
end

%% sweep the mode
for modei=1:length(ModeSelect)
    Residue = ResidueAll{modei};
    ZmVal = ZmValAll{modei};
    SigmaMag = abs(real(MdMode(ModeSelect(modei))))*2*pi; %MdMode is in the unite of Hz, so needs to be changed to rad.

    for j=1:length(CIMR2)
        k=CIMR2(j).device;
        if ApparatusType{k} ~= 100
            IMR = SigmaMag/(SimplusGT.Frobenius_norm_dq(Residue(k))*SimplusGT.Frobenius_norm_dq(ZmVal(k)));
            IMR_o=IMR;
        if IMR<0.01
            IMR = log10(0.01);
        else
            IMR = log10(IMR);
        end

        if IMR<CIMR2(j).value
            CIMR2(j).value=IMR;
            CIMR2(j).mode = MdMode(ModeSelect(modei));
            CIMR2(j).value_orig=IMR_o;
        end
        end
    end
end

% figure(3);
% clf;
% bar([CIMR(:).value]);
%Main_K_Plot(CIMR2,1);
title('Small-Signal System Strength Heatmap');