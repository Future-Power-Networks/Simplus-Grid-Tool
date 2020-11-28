function ExcelWrite(N_Bus,N_Device,DeviceType,DeviceStateStr,DeviceInputStr,DeviceOutputStr,ZbusStateStr, GminSS, GsysDSS)
%this function is to write states, devices, axes to select for
%further Greybox analysis.
% Author: Yue Zhu
filename = 'GreyBoxConfig.xlsx';

%%
%*********States sheet write
%Based on GsysDSS.
[~,D]=eig(GsysDSS.A,GsysDSS.E);
D=D/(2*pi);
Mode=diag(D);
ModeNum = length(Mode);

xlswrite(filename,{'Select states for state participation analysis. Only 1 mode should be selected.'},'State-PF','A1');
xlswrite(filename,{'write "1" for for selection, others for not'},'State-PF','A2');
StartSpace='A6';
StateSheet(1,1) = {'Device'};
StateSheet(1,2) = {'State'};
StateSheet(1,3) = {'Select'};
index = 2;
for k = 1:N_Bus
    if DeviceType{k} <= 89 %devices)
        DeviceName=strcat('Device',num2str(k));
        StateName = DeviceStateStr{k};
        StateNum = length(DeviceStateStr{k});
        StateSheet(index,1) = {DeviceName};
        for j = 1: StateNum
            StateSheet(index,2) = {StateName{j}};
            StateSheet(index,3) = {0};
            index = index+1;
        end
    else % floating bus, infinite bus...
    end
end
StateSheet(index,1) = {'Line'};
for i = 1: length(ZbusStateStr)
    StateSheet{index,2} = ZbusStateStr{i};
    StateSheet(index,3) = {0};
    index = index + 1;
end

StateSheet(1,5) = {'Mode'};
StateSheet(1,6) = {'Value'};
StateSheet(1,7) = {'Select'};
index = 2;
for i=1:ModeNum
    ModeName = strcat('Mode',num2str(i));
    StateSheet(index,5) = {ModeName};
    StateSheet(index,6) = {num2str(Mode(i),'%.2f')};
    StateSheet(index,7) = {0};
    index=index+1;
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

xlswrite(filename,{'Select devices and mode for bode-plot and Greybox Layer1&2 analysis'},'Impedance-PF','A1');
xlswrite(filename,{'write "1" for for selection, others for not'},'Impedance-PF','A2');

%*** Layer1&2 device select
ImpedanceSheet(1,1) = {'Device selection for Layer1&2'};
ImpedanceSheet(2,1) = {'Device'};
ImpedanceSheet(2,2) = {'Select'};
index=3;
for k = 1:N_Bus
        if DeviceType{k} <= 89 %devices)
            DeviceName=strcat('Device',num2str(k));
            ImpedanceSheet(index,1) = {DeviceName};
            ImpedanceSheet(index,2) = {0};
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
ImpedanceSheet(3,5) = {0};
ImpedanceSheet(4,5) = {0};
ImpedanceSheet(5,5) = {0};
ImpedanceSheet(6,5) = {0};

%*** Mode select.
ImpedanceSheet(1,7) = {'Mode selection for Layer 1&2&3'};
ImpedanceSheet(2,7) = {'Mode'};
ImpedanceSheet(2,8) = {'Value'};
ImpedanceSheet(2,9) = {'Select'};
index=3;
for i=1:ModeNum
    ModeName = strcat('Mode',num2str(i));
    ImpedanceSheet(index,7)={ModeName};
    ImpedanceSheet(index,8)={num2str(Mode(i),'%.2f')};
    ImpedanceSheet(index,9)={0};
    index=index+1;
end

%***Layer3 device Select
ImpedanceSheet(1,11)={'Device selection for Layer3'};
ImpedanceSheet(2,11)={'Device'};
ImpedanceSheet(2,12)={'Select'};
index=3;
for k = 1:N_Bus
        if DeviceType{k} <= 89 %devices)
            DeviceName=strcat('Device',num2str(k));
            ImpedanceSheet(index,11)= {DeviceName};
            ImpedanceSheet(index,12)= {0};
            index=index+1;
        else % floating bus, infinite bus...
        end
end
xlswrite(filename,ImpedanceSheet,'Impedance-PF',StartSpace);

%%
%%
end