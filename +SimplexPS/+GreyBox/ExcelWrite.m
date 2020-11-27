function ExcelWrite(N_Bus,N_Device,DeviceType,DeviceStateStr,DeviceInputStr,DeviceOutputStr,ZbusStateStr, GminSS)
%this function is to write states, devices, axes to select for
%further Greybox analysis.
% Author: Yue Zhu

%%
%States sheet write
fprintf('Writing states... \n');
filename = 'GreyBoxConfig.xlsx';
xlswrite(filename,{'Select states for state participation analysis. Only 1 mode should be selected.'},'States','A1');
xlswrite(filename,{'write "1" for for selection, others for not'},'States','A2');
index = 4;
for k = 1:N_Bus
    if DeviceType{k} <= 89 %devices)
        DeviceName=strcat('Device',num2str(k));
        StateName = DeviceStateStr{k};
        StateNum = length(DeviceStateStr{k});
        Space1 = strcat('A',num2str(index));
        Space2 = strcat('B',num2str(index));
        xlswrite(filename,{DeviceName},'States',Space1);
        xlswrite(filename,{'Select'},'States',Space2);
        for j = 1: StateNum
            Space1 = strcat('A',num2str(index+j));
            xlswrite(filename,{StateName{j}},'States',Space1);
        end
        index = index + 1 + StateNum;
    else % floating bus, infinite bus...
    end
end

%%
%Devices sheet write
fprintf('Writing devices for Layer1&2... \n');
xlswrite(filename,{'Select devices and mode for bode-plot and Greybox Layer1&2 analysis'},'Devices','A1');
xlswrite(filename,{'write "1" for for selection, others for not'},'Devices','A2');
xlswrite(filename,{'Devices'},'Devices','A4');
xlswrite(filename,{'Select'},'Devices','B4');
index=4;
for k = 1:N_Bus
        if DeviceType{k} <= 89 %devices)
            DeviceName=strcat('Device',num2str(k));
            Space = strcat('A',num2str(index+1));
            xlswrite(filename,{DeviceName},'Devices',Space);
            index=index+1;
        else % floating bus, infinite bus...
        end
end

xlswrite(filename,{'Bodeplot axis'},'Devices','D4');
xlswrite(filename,{'d-d'},'Devices','D5');
xlswrite(filename,{'d-q'},'Devices','D6');
xlswrite(filename,{'q-d'},'Devices','D7');
xlswrite(filename,{'q-q'},'Devices','D8');
xlswrite(filename,{'Select'},'Devices','E4');

xlswrite(filename,{'Mode'},'Devices','G4');
xlswrite(filename,{'Damping Coefficient'},'Devices','H4');
xlswrite(filename,{'Natural Frequency'},'Devices','I4');
xlswrite(filename,{'Select'},'Devices','J4');

[~,D]=eig(GminSS.A);
D=D/(2*pi);
Mode=diag(D);
ModeReal = real(Mode);
ModeImag = imag(Mode);
Num = length(Mode);
for i=1:Num
    ModeName = strcat('Mode',num2str(i));
    Space1 = strcat('G',num2str(i+4));
    Space2 = strcat('H',num2str(i+4));
    Space3 = strcat('I',num2str(i+4));
    xlswrite(filename,{ModeName},'Devices',Space1);
    xlswrite(filename,{num2str(ModeReal(i),'%.2f')},'Devices',Space2);
    xlswrite(filename,{num2str(ModeImag(i),'%.2f')},'Devices',Space3);
end


%%
%parameters sheet write
fprintf('Writing devices for Layer3... \n');
xlswrite(filename,{'Select devices for Layer3 analysis.'},'Parameters','A1');
xlswrite(filename,{'write "1" for for selection, others for not'},'Parameters','A2');
xlswrite(filename,{'Devices'},'Parameters','A4');
xlswrite(filename,{'Select'},'Parameters','B4');
index=4;
for k = 1:N_Bus
        if DeviceType{k} <= 89 %devices)
            DeviceName=strcat('Device',num2str(k));
            Space = strcat('A',num2str(index+1));
            xlswrite(filename,{DeviceName},'Parameters',Space);
            index=index+1;
        else % floating bus, infinite bus...
        end
end
%%
fprintf('GreyboxConfig.xlsx is now ready. Plese open the file and select the states and devices interested.\n');
fprintf('After selection, save the excel file and run GreyBoxAnalysis.m.\n');
%%
end