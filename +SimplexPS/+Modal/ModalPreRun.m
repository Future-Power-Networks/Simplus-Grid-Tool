% this .m file is used to prepare an excel file for user to config
% Participation / Sensitivity analysis. If the file existed, this program
% will clear the original contents. The excel file is built based on the
% running results of the toolbox.

% Author: Yue Zhu

FileModal=[cd '\NetworkConfigData\ModalConfig_' Name_Netlist '.xlsx'];%[cd '\ModalConfig.xlsx'];
%filename = 'ModalConfig.xlsx';
SimplexPS.Modal.ExcelPrep(FileModal); %create a new excel file, or clear old contents.
% write contents in the excel file.

% *For State PF analysis, Auto-Select will select all states at each device. It will
% also select two modes lower than 100Hz, with most small damping coefficient, but not
% around 0Hz (0.1Hz tollerence) or Fbase(1Hz tollerence).
% *For Impedance PF analysis, Auto-select will select all devices for
% Layer-1 analysis, d-d axis for bode plot, and device-1 for Layer-3
% analysis. As well, it will select two modes with lowest damping.
% Write 1 to enable auto select.
AutoSel = 0;

[AutoSelResult] = SimplexPS.Modal.ExcelWrite(N_Bus,N_Device,DeviceType,...
    DeviceStateStr,DeviceInputStr,DeviceOutputStr,ZbusStateStr, GminSS, GsysDSS, AutoSel, Fbase, FileModal);

fprintf('%s is now ready. Plese open the file and select the states and devices you are interested.\n',FileModal);
fprintf('After selection, save the excel file and run Modal Analysis.m.\n');

%winopen([Name_Netlist,'.xlsx']);
if AutoSel == 0
    winopen(FileModal);
end
if AutoSel == 1 && AutoSelResult == 1
    SimplexPS.Modal.ModalAnalysis
elseif AutoSel ==1 && AutoSelResult == 0
    error('Mode Auto-Selection failed. Please open ModalConfig.xlsx file to select the mode manually.')
else
end