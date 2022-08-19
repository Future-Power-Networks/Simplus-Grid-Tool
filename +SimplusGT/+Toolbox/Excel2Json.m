% This function converts an excel file to json file.
%
% Author(s): Rob Oldaker
%
% Modified by Yitong Li:
% # More comments are added.
% # The data operation of apparatuses is improved.
% # Partial address is added into file name.

function Excel2Json(FileName)

    fprintf('Convert "excel" to "json". \n')

    %
    % Get file name
    %
    RootPath = mfilename('fullpath');
    [RootPath,~,~]  = fileparts(RootPath);
    [RootPath,~,~]  = fileparts(RootPath);
    [RootPath,~,~]  = fileparts(RootPath);
    FilePath = fileparts(which(FileName));
    FolderName = erase(FilePath,RootPath);
    ExcelFile = [FolderName,'\',FileName];
    if strcmp(ExcelFile(1),'\')
        ExcelFile = ExcelFile(2:end);
    end
    
    
    %
    % Use a struct to hold the data
    %
    Data=struct;

    %
    % Basic settings
    %
    ExcelBasic = xlsread(ExcelFile,'Basic');
    Data.Basic = toBasic(ExcelBasic);

    %
    % Advanced settings
    %
    ExcelAdvance = xlsread(ExcelFile,'Advance');
    Data.Advance = toAdv(ExcelAdvance);

    %
    % Buses
    %
    ExcelBus = xlsread(ExcelFile,'Bus');
    Data.Bus=[];
    for i = 1:size(ExcelBus,1)
        Data.Bus = [Data.Bus toBus(ExcelBus(i,:))];
    end

    %
    % Network lines
    %
    ExcelNetworkLine = xlsread(ExcelFile,'NetworkLine');
    Data.NetworkLine=[];
    for i = 1:size(ExcelNetworkLine,1)
        Data.NetworkLine = [Data.NetworkLine toNetworkLine(ExcelNetworkLine(i,:))];
    end

    %
    % Network lines (IEEE)
    %
    ExcelNetworkLineIEEE = xlsread(ExcelFile,'NetworkLine_IEEE');
    Data.NetworkLineIEEE=[];
    for i = 4:size(ExcelNetworkLineIEEE,1)
        Data.NetworkLineIEEE = [Data.NetworkLineIEEE toNetworkLineIEEE(ExcelNetworkLineIEEE(i,:))];
    end

    % 
    % Apparatus
    %
    Wbase = Data.Basic.Fbase*2*pi;     % (rad/s), base angular frequency
    [ApparatusBus,ApparatusType,Para,NumApparatus] = SimplusGT.Toolbox.RearrangeListApparatus(ExcelFile,Wbase,ExcelBus);
    Data.Apparatus = [];
    for i = 1:NumApparatus
        Data.Apparatus = [Data.Apparatus toApparatus(ApparatusBus{i},ApparatusType{i},Para{i})];
    end
   
    % disp(Data);
    JsonFile = replace(ExcelFile,'.xlsx','.json');       % For xlsx case
    JsonFile = replace(JsonFile,'.xlsm','.json');           % For xlsm case
    if ~strcmp(JsonFile,ExcelFile)
        SimplusGT.SaveAsJsonToFile(Data,JsonFile);
        fprintf(['Successfully convert "', FileName, '" into ".json" type \n']);
    end

end

function [bus]=toBus(dBus)
    bus = struct;
    bus.BusNo   = dBus(1);
    bus.BusType = dBus(2);
    bus.Voltage = dBus(3);
    bus.Theta   = dBus(4);
    bus.PGi     = dBus(5);
    bus.QGi     = dBus(6);
    bus.PLi     = dBus(7);
    bus.QLi     = dBus(8);
    bus.Qmin    = dBus(9);
    bus.Qmax    = dBus(10);
    bus.AreaNo  = dBus(11);
    bus.AcDc    = dBus(12);
end

function [basic]=toBasic(b)
    basic = struct;
    basic.Fs = b(1);
    basic.Fbase = b(2);   % (Hz), base frequency
    basic.Sbase = b(3);   % (VA), base power
    basic.Vbase = b(4);   % (V), base voltage        
end

function [Advance]=toAdv(a)
    Advance = struct;
    Advance.DiscretizationMethod        = a(1);
    Advance.LinearizationTimes          = a(2);
    Advance.DiscretizationDampingFlag   = a(3);
    Advance.DirectFeedthrough           = a(4);
    Advance.PowerFlowAlgorithm    	    = a(5);
    Advance.EnableCreateSimulinkModel	= a(6);
    Advance.EnablePlotPole           	= a(7);
    Advance.EnablePlotAdmittance     	= a(8);
    Advance.EnablePrintOutput       	= a(9);
    Advance.EnableParticipation         = a(10);
end

function [Apparatus]=toApparatus(busNo,type,para)
    Apparatus = struct;
    Apparatus.BusNo  = busNo;
    Apparatus.Type   = type;
    Apparatus.Para = para;
end

function [NetworkLine]=toNetworkLine(nLine)
    NetworkLine = struct;
    NetworkLine.FromBus     = nLine(1);
    NetworkLine.ToBus       = nLine(2);
    NetworkLine.R           = nLine(3);
    NetworkLine.wL          = nLine(4);
    NetworkLine.wC          = nLine(5);
    NetworkLine.G           = nLine(6);
    NetworkLine.TurnsRatio  = nLine(7);
end

function [NetworkLine]=toNetworkLineIEEE(nLine)
    NetworkLine = struct;
    NetworkLine.FromBus     = nLine(1);
    NetworkLine.ToBus       = nLine(2);
    NetworkLine.R           = nLine(3);
    NetworkLine.X           = nLine(4);
    NetworkLine.B           = nLine(5);
    NetworkLine.G           = nLine(6);
    NetworkLine.TurnsRatio  = nLine(7);
end