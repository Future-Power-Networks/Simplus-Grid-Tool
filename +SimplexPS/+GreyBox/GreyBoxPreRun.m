% this .m file is used to prepare an excel file for user to config 
% participation analysis. If the file existed, this program will clear the
% original contents. The excel file is built based on the running results
% of the toolbox.
% Author: Yue Zhu

SimplexPS.GreyBox.ExcelPrep; %create a new excel file, or clear old contents.
% write contents in the excel file.
SimplexPS.GreyBox.ExcelWrite(N_Bus,N_Device,DeviceType,...
    DeviceStateStr,DeviceInputStr,DeviceOutputStr,ZbusStateStr, GminSS, GsysDSS);
fprintf('GreyboxConfig.xlsx is now ready. Plese open the file and select the states and devices you are interested.\n');
fprintf('After selection, save the excel file and run GreyBoxAnalysis.m.\n');

