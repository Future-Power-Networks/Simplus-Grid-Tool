function [AutoSelResult] = ExcelWrite(N_Bus,N_Apparatus,ApparatusType,ApparatusStateStr,ApparatusInputStr,...
    ApparatusOutputStr,ZbusStateStr, GminSS, GsysDSS, AutoSel, Fbase, filename)
%this function is to write states, apparatuses, axes to select for
%further Modal analysis.
% Author: Yue Zhu

AutoSelResult = 0;
%%
%*********States sheet write
%Based on GsysDSS.
[GsysSS, IndexSS] = SimplusGT.dss2ss(GsysDSS);
[~,D]=eig(GsysSS.A);
D=D/(2*pi);
Mode=diag(D);
ModeNum = length(Mode);

xlswrite(filename,{'Select states for state participation analysis. Only 1 mode should be selected.'},'State-PF','A1');
xlswrite(filename,{'write "1" for for selection, others for not'},'State-PF','A2');
StartSpace='A6';
StateSheet(1,1) = {'Apparatus'};
StateSheet(1,2) = {'State'};
StateSheet(1,3) = {'Select'};
index = 2;
StateCount = 0;
for k = 1:N_Bus
    if ApparatusType{k} <= 89 %apparatuses)    
        ApparatusName=strcat('Apparatus',num2str(k));
        StateName = ApparatusStateStr{k};
        StateNum = length(ApparatusStateStr{k});
        StateSheet(index,1) = {ApparatusName};
        for j = 1: StateNum
            StateCount = StateCount +1;
            if ismember(StateCount,IndexSS)
                StateSheet(index,2) = {StateName{j}};
                StateSheet(index,3) = {1};
%                 if (AutoSel==1 && j == 1) || (AutoSel==1 && j == 2) %select 'epsilon', and 'id' for pf analysis for each apparatus.
%                     StateSheet(index,3) = {1};
%                 else
%                     StateSheet(index,3) = {0};
%                 end
                index = index+1;
            end
        end
    else % floating bus, infinite bus...
    end
end

StateSheet(index,1) = {'Line'};
for i = 1: length(ZbusStateStr)
    StateCount = StateCount +1;
    if ismember(StateCount,IndexSS)
        StateSheet{index,2} = ZbusStateStr{i};
        StateSheet(index,3) = {0};
        index = index + 1;
    end
end

StateSheet(1,5) = {'Mode'};
StateSheet(1,6) = {'Value'};
StateSheet(1,7) = {'Select'};
index = 2;
SmodeSel1=abs(real(Mode(1)));
IndexSel1=2;
SmodeSel2=abs(real(Mode(2)));
IndexSel2=3;
for i=1:ModeNum
    ModeName = strcat('Mode',num2str(i));
    StateSheet(index,5) = {ModeName};
    StateSheet(index,6) = {num2str(Mode(i),'%.2f')};
    %auto select two modes
    if AutoSel ==1 && i >= 3
        if abs(imag(Mode(i))) < 100 % below 100Hz mode
            if abs(abs(imag(Mode(i))) - 0) >= 0.1 % not around 0.
                if abs(abs(imag(Mode(i))) - Fbase) >= 1 % not around Fbase
                    if abs(real(Mode(i))) < SmodeSel1 || SmodeSel1 == 0
                        SmodeSel2 = SmodeSel1;
                        IndexSel2 = IndexSel1;
                        SmodeSel1 = abs(real(Mode(i)));
                        IndexSel1 = index;
                    elseif abs(real(Mode(i)))> SmodeSel1 && abs(real(Mode(i)))< SmodeSel2 ...
                            || SmodeSel2 == 0 || abs(abs(imag(Mode(i))) - Fbase) <=1
                        SmodeSel2 = abs(real(Mode(i)));
                        IndexSel2 = index;
                    end
                end
            end
        end
    end
    StateSheet(index,7) = {0};
    index=index+1;
end
if IndexSel1~=2 && IndexSel2 ~= 3 % auto select success
    AutoSelResult = 1;
    StateSheet(IndexSel1,7) = {1};
    StateSheet(IndexSel2,7) = {1};
elseif AutoSel==1
    AutoSelResult = 0;
    error('Mode Auto-Selection failed. Please open ModalConfig.xlsx file to select the mode manually.')
else
end




xlswrite(filename,StateSheet,'State-PF',StartSpace);

%%
%Impedance-PF sheet write
%Based on GminSS
[~,D]=eig(GminSS.A);
D=D/(2*pi);
Mode=diag(D);
ModeReal = real(Mode);
ModeImag = imag(Mode);
ModeNum = length(Mode);

StartSpace='A6';

xlswrite(filename,{'Select apparatuses and mode for bode-plot and Modal Layer1&2 analysis'},'Impedance-PF','A1');
xlswrite(filename,{'write "1" for for selection, others for not'},'Impedance-PF','A2');

%*** Layer1&2 apparatus select
ImpedanceSheet(1,1) = {'Apparatus selection for Layer1&2'};
ImpedanceSheet(2,1) = {'Apparatus'};
ImpedanceSheet(2,2) = {'Select'};
index=3;
for k = 1:N_Bus
        if ApparatusType{k} <= 89 %apparatuses)
            ApparatusName=strcat('Apparatus',num2str(k));
            ImpedanceSheet(index,1) = {ApparatusName};
            ImpedanceSheet(index,2) = {1};
            index=index+1;
        else % floating bus, infinite bus...
        end
end

%*** bodplot d-q select
ImpedanceSheet(1,4) ={'Bodeplot axis selection'};
ImpedanceSheet(2,4) = {'Axis'};
ImpedanceSheet(2,5) = {'Select'};
ImpedanceSheet(3,4) = {'d-d'};
ImpedanceSheet(4,4) = {'d-q'};
ImpedanceSheet(5,4) = {'q-d'};
ImpedanceSheet(6,4) = {'q-q'};
if AutoSel == 1
    ImpedanceSheet(3,5) = {1};
else 
    ImpedanceSheet(3,5) = {0};
end
ImpedanceSheet(4,5) = {0};
ImpedanceSheet(5,5) = {0};
ImpedanceSheet(6,5) = {0};

%*** Mode select.
ImpedanceSheet(1,7) = {'Mode selection for Layer 1&2&3'};
ImpedanceSheet(2,7) = {'Mode'};
ImpedanceSheet(2,8) = {'Value'};
ImpedanceSheet(2,9) = {'Select'};
index=3;

SmodeSel1=abs(real(Mode(1)));
IndexSel1=3;
SmodeSel2=abs(real(Mode(2)));
IndexSel2=4;
for i=1:ModeNum
    ModeName = strcat('Mode',num2str(i));
    ImpedanceSheet(index,7)={ModeName};
    ImpedanceSheet(index,8)={num2str(Mode(i),'%.2f')};
    
    if AutoSel ==1 && i >= 3
        if abs(imag(Mode(i))) < 100 % below 100Hz mode
            if abs(abs(imag(Mode(i))) - 0) >= 0.1 % not around 0.
                if abs(abs(imag(Mode(i))) - Fbase) >= 1 % not around Fbase
                    if abs(real(Mode(i))) < SmodeSel1 || SmodeSel1 == 0 || abs(abs(imag(Mode(i))) - Fbase) <=1
                        SmodeSel2 = SmodeSel1;
                        IndexSel2 = IndexSel1;
                        SmodeSel1 = abs(real(Mode(i)));
                        IndexSel1 = index;
                    elseif abs(real(Mode(i)))> SmodeSel1 && abs(real(Mode(i)))< SmodeSel2 ...
                            || SmodeSel2 == 0 || abs(abs(imag(Mode(i))) - Fbase) <=1
                        SmodeSel2 = abs(real(Mode(i)));
                        IndexSel2 = index;
                    end
                end
            end
        end
    end
    ImpedanceSheet(index,9)={0};
    index=index+1;
end
if AutoSel==1 && IndexSel1~=2 && IndexSel2 ~= 3 % auto select success
    AutoSelResult = 1;
    ImpedanceSheet(IndexSel1,9) = {1};
    ImpedanceSheet(IndexSel2,9) = {1};
elseif AutoSel==1
    AutoSelResult = 0;
    error('Mode Auto-Selection failed. Please open ModalConfig.xlsx file to select the mode manually.')
else
end


%***Layer3 apparatus Select
ImpedanceSheet(1,11)={'Apparatus selection for Layer3'};
ImpedanceSheet(2,11)={'Apparatus'};
ImpedanceSheet(2,12)={'Select'};
index=3;
IndexSel=1;
for k = 1:N_Bus
        if ApparatusType{k} <= 89 %apparatuses)
            ApparatusName=strcat('Apparatus',num2str(k));
            ImpedanceSheet(index,11)= {ApparatusName};
            if AutoSel == 1 && IndexSel==1
                ImpedanceSheet(index,12) = {1};
                IndexSel = IndexSel+1;
            else
                ImpedanceSheet(index,12) = {0};
            end
            index=index+1;
        else % floating bus, infinite bus...
        end
end
xlswrite(filename,ImpedanceSheet,'Impedance-PF',StartSpace);

%% Enabling sheet writing
EnableSheet(1,1) = {'Select the function you wish to enable'};
EnableSheet(2,1) = {'write "1" for for selection, others for not'};

EnableSheet(6,1) = {'function'};
EnableSheet(6,2) = {'Select'};
EnableSheet(7,1) = {'State-PF'};
EnableSheet(8,1) = {'Bodeplot'};
EnableSheet(9,1) = {'Impedance-PF Layer 1&2'};
EnableSheet(10,1) = {'Impedance-PF Layer 3'};
for i=7:10
    EnableSheet(i,2) = {1};
end
xlswrite(filename,EnableSheet,'Enabling','A1');


%%
%%
end