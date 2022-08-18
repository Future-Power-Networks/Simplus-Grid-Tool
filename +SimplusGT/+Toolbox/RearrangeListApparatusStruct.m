% This function re-arranges the netlist data of apparatuses.

% Author(s): Yitong Li, Yunjie Gu
% Modified by Rob

%% Notes
%
% The apparatus model is in load convention.
%
% The apparautus bus index, type index, and parameters are saved in three
% different cells. The reason of doing this is that, the number of
% apparatuses might not be same to the number of buses, when the power
% system has interlink apparatuses. In this case, the interlink apparatus
% will have two bus indice.
%
% "app" means apparatus.

function [AppBusCell,AppTypeCell,ParaCell,N_App] = RearrangeListApparatusStruct(InputData,W0,ListBus)

ListApp=[];
index=1;
ListAppBusChar={};
for i=1:size(InputData.Apparatus,1)
    aData = InputData.Apparatus(i);
    [ListAppArray,ListAppCharArray] = toListAppArray(aData);
    ListApp(i,:) = ListAppArray;
    ListAppBusChar{i} = ListAppCharArray;
    index=index+1;
end

%% Rearrange data
% Notes:
% This section will convert the ListAppBus and ListAppType to AppBusCell
% and AppTypeCell, i.e., array type to cell type.

[N_App,ColumnMax_Apparatus] = size(ListApp);
ListAppBus = ListApp(:,1);
ListAppType = ListApp(:,2);

% Get the apparatus bus in cell form
for n = 1:N_App
    if ~isnan(ListAppBus(n))  % If not NaN, then the bus is a scalar rather han an array in char type
        AppBusCell{n} = ListAppBus(n);
    else
        AppBusCell{n} = str2num(ListAppBusChar{n});
        [~,~,AreaType]= SimplusGT.Toolbox.CheckBus(AppBusCell{n}(1),ListBus);
        
        % Notes:
        % If the first bus is dc bus, then swap. This ensures the first one
        % is ac bus.
        if AreaType == 2 
            [AppBusCell{n}(1),AppBusCell{n}(2)] = deal(AppBusCell{n}(2),AppBusCell{n}(1));
        end
    end
end

% Get the apparatus type in cell form
for n = 1:N_App
    AppTypeCell{n} = ListApp(n,2);
end

% Add floating bus
N_Bus = length(ListBus(:,1));
for m = 1:N_Bus
    BusIndex = ListBus(m,1);
    ExistApparatus = SimplusGT.CellFind(AppBusCell,BusIndex);
    if isempty(ExistApparatus)
        % The bus has no apparatus, i.e., an ac or dc floating bus.
        N_App = N_App+1;
        AppBusCell{N_App} = BusIndex;   % Add a new bus index
        [~,~,AreaType] = SimplusGT.Toolbox.CheckBus(BusIndex,ListBus);  % Check if an ac or dc bus
        if AreaType == 1
            AppTypeCell{N_App} = 100;       % Ac floating bus
        elseif AreaType == 2
            AppTypeCell{N_App} = 1100;      % Dc floating bus
        else
            error(['Error: Error AreaType.']);
        end
    else
        % The bus has an apparatus already, no need to do anything then.
    end
end

% Error check
if (ColumnMax_Apparatus>13)
    error(['Error: Apparatus data overflow.']); 
end

[~,ModeBus] = SimplusGT.CellMode(AppBusCell);
if ModeBus~=1
    error(['Error: For each bus, one and only one apparatus has to be connected.']);
end

%% Default AC apparatus data
% ======================================
% Synchronous generator
% ======================================
Para0000.J  = 3.5;
Para0000.D  = 1;
Para0000.wL = 0.05;
Para0000.R  = 0.01;
Para0000.w0 = W0;

% ======================================
% Grid-following VSI (PLL-controlled)
% ======================================
% Dc link
Para0010.V_dc       = 2.5;
Para0010.C_dc       = 1.25;         % 2*0.1*Para0010.V_dc^2;
Para0010.f_v_dc     = 5;            % (Hz) bandwidth, vdc

% Ac filter
Para0010.wLf        = 0.03;
Para0010.R          = 0.01;

% PLL
Para0010.f_pll      = 5;            % (Hz) bandwidth, PLL
Para0010.f_tau_pll  = 300;          % (Hz) bandwidth, PLL low pass filter

% Current loop
Para0010.f_i_dq     = 600;      	% (Hz) bandwidth, idq
Para0010.w0         = W0;   

% ======================================
% Grid-forming VSI (Droop-Controlled)
% ======================================
Para0020.wLf    =0.05;
Para0020.Rf     =0.05/5;
Para0020.wCf    =0.02;
Para0020.wLc    =0.01;
Para0020.Rc     =0.01/5;
Para0020.Xov    =0.01;
Para0020.Dw     =0.05;
Para0020.fdroop =5;    % (Hz) droop control bandwidth
Para0020.fvdq   =300;   % (Hz) vdc bandwidth
Para0020.fidq   =600;   % current control bandwidth
Para0020.w0     = W0;

% ======================================
% Ac infinite bus (short-circuit in small-signal)
% ======================================
Para0090 = [];

% ======================================
% Ac floating bus (open-circuit)
% ======================================
Para0100 = [];

%% Default DC apparatus data
% ======================================
% Grid-feeding buck
% ======================================
Para1010.Vdc  = 2;
Para1010.Cdc  = 0.8;
Para1010.wL   = 0.05;
Para1010.R    = 0.05/5;
Para1010.fi   = 600;
Para1010.fvdc = 5;
Para1010.w0   = W0;

% ======================================
% Dc infinite bus (short-circuit in small-signal)
% ======================================
Para1090 = [];

% ======================================
% Dc floating bus (open-circuit)
% ======================================
Para1100 = [];

%% Default hybrid apparatus data
% ======================================
% Interlink ac-dc converter
% ======================================
Para2000.C_dc   = 1.6;
Para2000.wL_ac  = 0.05;
Para2000.R_ac   = 0.01;
Para2000.wL_dc  = 0.02;
Para2000.R_dc   = 0.02/5;
Para2000.fidq   = 600;
Para2000.fvdc   = 5;
Para2000.fpll   = 5;
Para2000.w0     = W0;   

%% Re-arrange apparatus data

% Find the index of user-defined data
ListApparatus_NaN = isnan(ListApp(:,3:ColumnMax_Apparatus));    % Find NaN
[row,column] = find(ListApparatus_NaN == 0);     
column = column+2;

% Initialize the apparatus parameters by default parameters
for i = 1:N_App
    AppBus   = AppBusCell{i};
    AppType  = AppTypeCell{i};
    switch floor(AppType/10)
        % ### AC apparatuses
        case 0     
            ParaCell{i} = Para0000;     % Synchronous machine
        case 1
            ParaCell{i} = Para0010;     % Grid-following inverter
      	case 2
            ParaCell{i} = Para0020;     % Grid-forming inverter
        case 3
            % Yue's Full-Order Machine
        case 9
            ParaCell{i} = Para0090;     % Ac inifnite bus
        case 10
            ParaCell{i} = Para0100;     % Ac floating bus, i.e., no apparatus
        
        % ### DC apparatuses
        case 101
            ParaCell{i} = Para1010;     % Grid-following buck
        case 109
            ParaCell{i} = Para1090;     % Dc infinite bus
        case 110
            ParaCell{i} = Para1100;     % Ac floating bus, i.e., no apparatus
            
        % ### Hybrid ac-dc apparatuses
        case 200
            ParaCell{i} = Para2000;     % Interlinking ac-dc converter
            
        % ### Error check
        otherwise
            error(['Error: apparatus type, bus ' num2str(AppBus) ' type ' num2str(AppType) '.']);
    end
end

%% The re-order can only be done here
% For a single-area pure ac system, we re-order the apparatuses 
AreaTypeCheck = find(ListBus(:,12) == 2, 1);
AreaNoCheck = find(ListBus(:,11) == 2, 1);
if isempty(AreaTypeCheck) && isempty(AreaNoCheck)
    for j = 1:N_App
        [~,AppIndex] = SimplusGT.CellFind(AppBusCell,j);
        AppBusCellNew{j}    = AppBusCell{AppIndex};
        AppTypeCellNew{j}   = AppTypeCell{AppIndex};
        ParaCellNew{j}      = ParaCell{AppIndex};
    end
    AppBusCell  = AppBusCellNew;
    AppTypeCell = AppTypeCellNew;
    ParaCell    = ParaCellNew;
end

end

%%
function [a,BusNoStr] = toListAppArray(aData)
a=nan(1,13);
BusNoStr = '';
if length(aData.BusNo)==1
    a(1) = aData.BusNo(1);
else
    a(1) = NaN;
    BusNoStr = '';
    for idx = 1:length(aData.BusNo)
        BusNoStr=strcat(BusNoStr,sprintf('%d',aData.BusNo(idx)));
        if idx<length(aData.BusNo)
            BusNoStr = strcat(BusNoStr,',');
        end
    end
end
a(2) = aData.Type;
if floor(aData.Type/10) == 0
    a(3) = aData.Para.J;
    a(4) = aData.Para.D;
    a(5) = aData.Para.wL;
    a(6) = aData.Para.R;
    a(7) = NaN;
    a(8) = NaN;
    a(9) = NaN;        
elseif floor(aData.Type/10) == 1
    a(3) = aData.Para.V_dc;
    a(4) = aData.Para.C_dc;
    a(5) = aData.Para.wLf;
    a(7) = aData.Para.f_v_dc;
    a(6) = aData.Para.R;
    a(8) = aData.Para.f_pll;
    a(9) = aData.Para.f_i_dq;
elseif floor(aData.Type/10) == 2
    a(3) = aData.Para.wLf;
    a(4) = aData.Para.Rf;
    a(5) = aData.Para.wCf;
    a(6) = aData.Para.wLc;
    a(7) = aData.Para.Rc;
    a(8) = aData.Para.Xov;
    a(9) = aData.Para.wCf;
    a(10) = aData.Para.Dw;
    a(11) = aData.Para.fdroop;
    a(12) = aData.Para.fvdq;
    a(13) = aData.Para.fidq;
elseif floor(aData.Type/10)==101
    a(3) = aData.Para.Vdc;
    a(4) = aData.Para.Cdc;
 	a(5) = aData.Para.wL;
    a(6) = aData.Para.R;
    a(7) = aData.Para.fi;
    a(8) = aData.Para.fvdc;
elseif floor(aData.Type/10)==200
    a(3) = aData.Para.C_dc;
    a(4) = aData.Para.wL_ac;
    a(5) = aData.Para.R_ac;
    a(6) = aData.Para.wL_dc;
    a(7) = aData.Para.R_dc;
    a(8) = aData.Para.fidq;
    a(9) = aData.Para.fvdc;
    a(10)= aData.Para.fpll;
end
end