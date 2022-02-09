file = 'UserData.xlsm';
convertXlsmToJson(file);
file = 'Examples\HybridPowerSystem\Hybrid_test_v1.xlsx';
convertXlsmToJson(file);
file = 'Examples\HybridPowerSystem\Hybrid_test_v2.xlsx';
convertXlsmToJson(file);
file = 'Examples\DcPowerSystem\SingleMachineInfiniteBus\GfdBuckInfiniteBus.xlsx';
convertXlsmToJson(file);
file = 'Examples\AcPowerSystem\SingleApparatusInfiniteBus\SgInfiniteBus.xlsx';
convertXlsmToJson(file);
file = 'Examples\AcPowerSystem\SingleApparatusInfiniteBus\GflInverterInfiniteBus.xlsx';
convertXlsmToJson(file);
file = 'Examples\AcPowerSystem\SingleApparatusInfiniteBus\GfmInverterInfiniteBus.xlsx';
convertXlsmToJson(file);
file = 'Examples\AcPowerSystem\NET_NYPS_68Bus\NETS_NYPS_68Bus.xlsx';
convertXlsmToJson(file);
file = 'Examples\AcPowerSystem\IEEE_57Bus\IEEE_57Bus.xlsx';
convertXlsmToJson(file);
file = 'Examples\AcPowerSystem\IEEE_30Bus\IEEE_30Bus.xlsx';
convertXlsmToJson(file);
file = 'Examples\AcPowerSystem\IEEE_14Bus\IEEE_14Bus.xlsx';
convertXlsmToJson(file);

function convertXlsmToJson(file)

    %
    % Use a struct to hold the data
    %
    data=struct;

    %
    % Basic settings
    %
    dataBasic = xlsread(file,'Basic');
    data.basic = toBasic(dataBasic);

    %
    % Advanced settings
    %
    dataAdv = xlsread(file,'Advance');
    data.adv = toAdv(dataAdv);

    %
    % Buses
    %
    [ListBus,N_Bus] = SimplusGT.Toolbox.RearrangeListBus(file);
    data.buses=[];
    for i = 1:N_Bus
        data.buses = [data.buses toBus(ListBus(i,:))];
    end

    %
    % Network lines
    %
    ListLine = xlsread(file,'NetworkLine');
    data.networkLines=[];
    sx = size(ListLine);
    for i = 1:size(ListLine,1)
        data.networkLines = [data.networkLines toNetworkLine(ListLine(i,:))];
    end

    %
    % Network lines (IEEE)
    %
    ListLine = xlsread(file,'NetworkLine_IEEE');
    data.networkLinesIEEE=[];
    for i = 4:size(ListLine,1)
        data.networkLinesIEEE = [data.networkLinesIEEE toNetworkLineIEEE(ListLine(i,:))];
    end

    % 
    % Apparatus
    %
    Wbase = data.basic.Fbase*2*pi;     % (rad/s), base angular frequency
    [ApparatusBus,ApparatusType,Para,N_Apparatus] = SimplusGT.Toolbox.RearrangeListApparatus(file,Wbase,ListBus);    
    data.apparatus = [];
    for i = 1:N_Apparatus
        data.apparatus = [data.apparatus toApparatus(ApparatusBus{i},ApparatusType{i},Para{i})];
    end
   
    disp(data);
    json = jsonencode(data,'PrettyPrint',true,'ConvertInfAndNaN',false);
    % NaN and Inf are not valid json literals so put them in quotes to make it valid json
    json = strrep(json,'-Inf','"-Inf"');
    json = strrep(json,'NaN','"NaN"');
    json = strrep(json,'Inf','"Inf"');

    jsonFile = replace(file,'.xlsx','.json');
    jsonFile = replace(jsonFile,'.xlsm','.json');
    if ~strcmp(jsonFile,file)
        fid = fopen(jsonFile,'wt');
        fprintf(fid, json, + '\n');
        fclose(fid);
        fprintf('Successfully converted file %s\n',file);
    end

    function [bus]=toBus(dBus)
        bus = struct;
        bus.busNo   = dBus(1)
        bus.busType = dBus(2);
        bus.VSp     = dBus(3);
        bus.theta   = dBus(4);
        bus.PGi     = dBus(5);
        bus.QGi     = dBus(6);
        bus.PLi     = dBus(7);
        bus.QLi     = dBus(8);
        bus.Qmin    = dBus(9);
        bus.Qmax    = dBus(10);
        bus.areaNo  = dBus(11);
        bus.ACDC    = dBus(12);
    end

    function [basic]=toBasic(b)
        basic = struct;
        basic.Fs = b(1);
        basic.Fbase = b(2);   % (Hz), base frequency
        basic.Sbase = b(3);   % (VA), base power
        basic.Vbase = b(4);   % (V), base voltage        
    end

    function [adv]=toAdv(a)
        adv = struct;
        adv.discretizationMethod        = a(1);
        adv.linearizationTimes          = a(2);
        adv.discretizationDampingFlag   = a(3);
        adv.directFeedThrough           = a(4);
        adv.powerFlowAlgorithm    	    = a(5);
        adv.enableCreateSimulinkModel	= a(6);
        adv.enablePlotPole           	= a(7);
        adv.enablePlotAdmittance     	= a(8);
        adv.enablePrintOutput       	= a(9);
        adv.enableParticipation         = a(10);
    end

    function [apparatus]=toApparatus(busNo,type,params)
        apparatus = struct;
        apparatus.busNo  = busNo;
        apparatus.type   = type;
        apparatus.params = params;
    end

    function [networkLine]=toNetworkLine(nLine)
        networkLine = struct;
        networkLine.fromBus     = nLine(1);
        networkLine.toBus       = nLine(2);
        networkLine.R           = nLine(3);
        networkLine.wL          = nLine(4);
        networkLine.wC          = nLine(5);
        networkLine.G           = nLine(6);
        networkLine.turnsRatio  = nLine(7);
    end

    function [networkLine]=toNetworkLineIEEE(nLine)
        networkLine = struct;
        networkLine.fromBus     = nLine(1);
        networkLine.toBus       = nLine(2);
        networkLine.R           = nLine(3);
        networkLine.X           = nLine(4);
        networkLine.B           = nLine(5);
        networkLine.G           = nLine(6);
        networkLine.turnsRatio  = nLine(7);
    end
end