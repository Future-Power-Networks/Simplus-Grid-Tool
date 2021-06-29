% This .m file is used to prepare an excel file for user to configure
% participation/sensitivity analysis. If the file existed, this program
% will clear the original contents. The excel file is built based on the
% running results of the toolbox.

% Author: Yue Zhu

% Warning added by Yitong:
% Running this function needs to close all excel forms. Then, the
% UserData.xlsx will be re-opened. This might influence the unsaved work of
% users.

% Generated variables: 
% UserData Pure: pure name of the user data excel file
% UserdataModal: full path of the modal config excel file
%%
% Set the file saving location and change the suffix
% this way avoids using  'cd', in case 'cd' is not set correctly at the
% toolbox.
UserdataModal=mfilename('fullpath');
pos_v =  strfind(UserdataModal,filesep);
UserdataModal = UserdataModal(1:pos_v(end-2));
clear pos_v;

UserDataPure = UserData;
UserDataPure = strrep(UserDataPure,'.xlsm','');
UserDataPure = strrep(UserDataPure,'.xlsx','');
UserDataPure = strrep(UserDataPure,'.xls','');
UserdataModal = [UserdataModal 'Examples\ParticipationAnalysis\','ModalConfig_', UserDataPure,'.xlsx'];

%%
% Set file name
%FileModal=[cd '\Examples\ParticipationAnalysis\', UserData_Modal];
SimplusGT.Modal.ExcelPrep(UserdataModal); %create a new excel file, or clear old contents.

%%
% write contents in the excel file.

% *For State PF analysis, Auto-Select will select all states at each apparatus. It will
% also select two modes lower than 100Hz, with most small damping coefficient, but not
% around 0Hz (0.1Hz tollerence) or Fbase(1Hz tollerence).
% *For Impedance PF analysis, Auto-select will select all apparatuses for
% Layer-1 analysis, d-d axis for bode plot, and apparatus-1 for Layer-3
% analysis. As well, it will select two modes with lowest damping.
% Write 1 to enable auto select.
MdAutoSel = 1;

[AutoSelResult] = SimplusGT.Modal.ExcelWrite(N_Apparatus,ApparatusType,ApparatusBus,...
    ApparatusStateStr,ZbusStateStr, GminSS, GsysDSS, MdAutoSel, Fbase, UserdataModal);

fprintf('%s is now ready.\nPlease open the file and select the states and apparatuses you are interested.\n',UserdataModal);
fprintf('After selection, save the excel file and run SimplusGT.Modal.ModalAnalysis.m.\n');
winopen(UserData);
if MdAutoSel == 0
    winopen(UserdataModal);
end
if MdAutoSel == 1 && AutoSelResult == 1
    SimplusGT.Modal.ModalAnalysis % run modal analysis
elseif AutoSelResult == 0  % do nothing
elseif MdAutoSel ==1 && AutoSelResult == 0
    error(['Error: Mode Auto-Selection failed. Please open ModalConfig.xlsx file to select the mode manually.'])
else
end

%winopen(UserdataModal);