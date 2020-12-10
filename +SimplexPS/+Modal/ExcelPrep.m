function ExcelPrep(FileModal)
%this function is to create a new excel file for users to config Modal
%analysis. If the ModalCofig.xlsx exists, this funciton will wipe it up,
%and rewrtie all the contents.
%further Modal analysis.
% Sheet modification is based on Microsoft Automation Server.
% Author: Yue Zhu

%fprintf('Preparing Excel file for user configuration...\n');
%pwd - current folder location.
%delete(FileModal); %delete the old file.
new=0;
try
    Excel = actxGetRunningServer('Excel.Application');
catch
    Excel = actxserver('Excel.Application');
end
if exist(FileModal,'file')
else
    eWorkbook = Excel.Workbooks.Add; %add a workbook
    SaveAs(eWorkbook,FileModal);
    Close(eWorkbook);
    new=1;
end
[~, sheet_names] = xlsfinfo(FileModal);
set(Excel, 'Visible', 0); %application invisible
set(Excel,'DisplayAlerts',0); %mute all alerts
Workbooks = Excel.Workbooks;  % Get a handle to Excel's Workbooks
Workbook=Workbooks.Open(FileModal); % Open an Excel Workbook and activate it
Sheets = Excel.ActiveWorkBook.Sheets;
if new ==1 
    Sheets.Add([], Sheets.Item(Sheets.Count));
end
SheetCount = Sheets.Count;
index_adjust = 0;
Sheets.Add([], Sheets.Item(SheetCount)); %add a new sheet after the last sheet.
Sheets.Add([], Sheets.Item(SheetCount+1));
Sheets.Add([], Sheets.Item(SheetCount+2));
for i=1:SheetCount %delete old sheets
    current_sheet = get(Sheets, 'Item', (i-index_adjust));
    invoke(current_sheet, 'Delete')
    %fprintf('Worksheet called %s deleted. \n',sheet_names{i});
    index_adjust = index_adjust +1;
end

%rename new blank sheets
Sheets.Item(1).Name = 'State-PF';
%fprintf('Worksheet called %s added. \n',sheet_names{1});
Sheets.Item(2).Name = 'Impedance-PF';
%fprintf('Worksheet called %s added. \n',sheet_names{2});
Sheets.Item(3).Name = 'Enabling';
%fprintf('Worksheet called %s added. \n',sheet_names{3});

Workbook.Save;
Workbooks.Close;
invoke(Excel, 'Quit');
delete(Excel);
end
