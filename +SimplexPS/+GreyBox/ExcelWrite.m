function ExcelWrite(N_Bus,N_Device,DeviceType,DeviceStateStr,DeviceInputStr,DeviceOutputStr,ZbusStateStr, GminSS)
%this function is to write states, devices, axes to select for
%further Greybox analysis.
% Author: Yue Zhu

%%
%States sheet write
fprintf('Writing values in state-PF sheet... \n');
filename = 'GreyBoxConfig.xlsx';
xlswrite(filename,{'Select states for state participation analysis. Only 1 mode should be selected.'},'State-PF','A1');
xlswrite(filename,{'write "1" for for selection, others for not'},'State-PF','A2');
IndexStart=6;

index = IndexStart;
for k = 1:N_Bus
    if DeviceType{k} <= 89 %devices)
        DeviceName=strcat('Device',num2str(k));
        StateName = DeviceStateStr{k};
        StateNum = length(DeviceStateStr{k});
        Space1 = strcat('A',num2str(index));
        Space2 = strcat('B',num2str(index));
        xlswrite(filename,{DeviceName},'State-PF',Space1);
        xlswrite(filename,{'Select'},'State-PF',Space2);
        for j = 1: StateNum
            Space1 = strcat('A',num2str(index+j));
            xlswrite(filename,{StateName{j}},'State-PF',Space1);
        end
        index = index + 1 + StateNum;
    else % floating bus, infinite bus...
    end
end

%%
%Impedance-PF sheet write
fprintf('Writing values in impdance-PF sheet... \n');
xlswrite(filename,{'Select devices and mode for bode-plot and Greybox Layer1&2 analysis'},'Impedance-PF','A1');
xlswrite(filename,{'write "1" for for selection, others for not'},'Impedance-PF','A2');

%***
xlswrite(filename,{'Device selection for Layer1&2'},'Impedance-PF',['A',num2str(IndexStart-1)]);
xlswrite(filename,{'Device'},'Impedance-PF',['A',num2str(IndexStart)]);
xlswrite(filename,{'Select'},'Impedance-PF',['B',num2str(IndexStart)]);
index=IndexStart;
for k = 1:N_Bus
        if DeviceType{k} <= 89 %devices)
            DeviceName=strcat('Device',num2str(k));
            xlswrite(filename,{DeviceName},'Impedance-PF',['A',num2str(index+1)]);
            xlswrite(filename,0,'Impedance-PF',['B',num2str(index+1)]);
            index=index+1;
        else % floating bus, infinite bus...
        end
end

%***
index=IndexStart;
xlswrite(filename,{'Bodeplot axis'},'Impedance-PF',['D',num2str(index)]);
xlswrite(filename,{'d-d'},'Impedance-PF',['D',num2str(index+1)]);
xlswrite(filename,{'d-q'},'Impedance-PF',['D',num2str(index+2)]);
xlswrite(filename,{'q-d'},'Impedance-PF',['D',num2str(index+3)]);
xlswrite(filename,{'q-q'},'Impedance-PF',['D',num2str(index+4)]);
xlswrite(filename,{'Select'},'Impedance-PF',['E',num2str(index)]);
xlswrite(filename,0,'Impedance-PF',['E',num2str(index+1)]);
xlswrite(filename,0,'Impedance-PF',['E',num2str(index+2)]);
xlswrite(filename,0,'Impedance-PF',['E',num2str(index+3)]);
xlswrite(filename,0,'Impedance-PF',['E',num2str(index+4)]);

%***
xlswrite(filename,{'Mode selection for Layer 1&2&3'},'Impedance-PF',['G',num2str(index-1)]);
xlswrite(filename,{'Mode'},'Impedance-PF',['G',num2str(index)]);
xlswrite(filename,{'Complex value'},'Impedance-PF',['H',num2str(index)]);
%xlswrite(filename,{'Natural Frequency'},'Impedance-PF',['I',num2str(index)]);
xlswrite(filename,{'Select'},'Impedance-PF',['I',num2str(index)]);

[~,D]=eig(GminSS.A);
D=D/(2*pi);
Mode=diag(D);
ModeReal = real(Mode);
ModeImag = imag(Mode);
Num = length(Mode);

index=IndexStart;
for i=1:Num
    ModeName = strcat('Mode',num2str(i));
    xlswrite(filename,{ModeName},'Impedance-PF',['G',num2str(index+1)]);
    xlswrite(filename,{num2str(Mode(i),'%.2f')},'Impedance-PF',['H',num2str(index+1)]);
    xlswrite(filename,0,'Impedance-PF',['I',num2str(index+1)]);
    %xlswrite(filename,{num2str(ModeImag(i),'%.2f')},'Impedance-PF',['I',num2str(index+1)]);
    index=index+1;
end

%***
xlswrite(filename,{'Device selection for Layer3'},'Impedance-PF',['K',num2str(IndexStart-1)]);
xlswrite(filename,{'Device'},'Impedance-PF',['K',num2str(IndexStart)]);
xlswrite(filename,{'Select'},'Impedance-PF',['L',num2str(IndexStart)]);
index=IndexStart;
for k = 1:N_Bus
        if DeviceType{k} <= 89 %devices)
            DeviceName=strcat('Device',num2str(k));
            xlswrite(filename,{DeviceName},'Impedance-PF',['K',num2str(index+1)]);
            xlswrite(filename,0,'Impedance-PF',['L',num2str(index+1)]);
            index=index+1;
        else % floating bus, infinite bus...
        end
end

% %%
% %parameters sheet write
% fprintf('Writing devices for Layer3... \n');
% xlswrite(filename,{'Select devices for Layer3 analysis.'},'Parameters','A1');
% xlswrite(filename,{'write "1" for for selection, others for not'},'Parameters','A2');
% xlswrite(filename,{'Impedance-PF'},'Parameters','A4');
% xlswrite(filename,{'Select'},'Parameters','B4');
% index=4;
% for k = 1:N_Bus
%         if DeviceType{k} <= 89 %devices)
%             DeviceName=strcat('Device',num2str(k));
%             Space = strcat('A',num2str(index+1));
%             xlswrite(filename,{DeviceName},'Parameters',Space);
%             index=index+1;
%         else % floating bus, infinite bus...
%         end
% end
%%
fprintf('GreyboxConfig.xlsx is now ready. Plese open the file and select the states and devices interested.\n');
fprintf('After selection, save the excel file and run GreyBoxAnalysis.m.\n');
%%
end