% This function re-arranges the netlist data of apparatuses.

% Author(s): Yitong Li, Yunjie Gu

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

function [AppBusCell,AppTypeCell,ParaCell,N_App] = RearrangeListApparatus(UserData,W0,ListBus)

%% New excel in xlsm format
% Check if the xlsm format is used
if strcmp(UserData,'UserData.xlsm')
    NewExcel = 1;
else
    NewExcel = 0;
end

%% Load data
[ListApp,ListAppChar]	 = xlsread(UserData,'Apparatus');

%% Prepare
if NewExcel == 1
    ListApp = ListApp(3:end,:);       % Remove the first two lines
    [rmax,cmax] = size(ListApp);
    ListAppNew = NaN(rmax/2,cmax);
    for r = 1:(rmax/2)
        ListAppNew(r,1) = ListApp(2*r-1,1);         % Move bus number
        ListAppNew(r,2) = ListApp(2*r-1,2);         % Move bus number
        ListAppNew(r,3:end) = ListApp(2*r,3:end);   % Move others
    end
   ListApp = ListAppNew;            % Update
end
ListAppBusChar = ListAppChar(:,1);
for k1 = 1:length(ListAppBusChar)
    if strcmpi(ListAppBusChar{k1},'Bus No.')
        break;
    end
end
ListAppBusChar = ListAppBusChar(k1+1:end);
if NewExcel == 1
    ListAppBusChar = ListAppBusChar(1:2:end,:);
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
if (ColumnMax_Apparatus>19)
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
% Battery Energy Storage System (Droop-Controlled)
% ======================================
Para0030.wLf      =0.05;
Para0030.Rf       =0.05/5;
Para0030.wCf      =0.02;
Para0030.wLc      =0.01;
Para0030.Rc       =0.01/5;
Para0030.Xov      =0.01;
Para0030.Dw       =0.05;
Para0030.fdroop   =5;    % (Hz) droop control bandwidth
Para0030.fvdq     =300;   % (Hz) vdc bandwidth
Para0030.fidq     =600;   % current control bandwidth
Para0030.w0       = W0;
Para0030.L_dc     =0.000125;
Para0030.C_dc     =0.224;
Para0030.v_ocv    =0.625; % open circuit voltage of the battery
Para0030.R_bat    =0.0625;
Para0030.v_dc_ref =1;
Para0030.fvdc     =50; % (Hz) vdc bandwidth
Para0030.fibat    =500; % (Hz) ibat bandwidth

% ======================================
% Photovoltaic (VFM-FFL)
% ======================================
Para0040.wLf      =0.05;
Para0040.Rf       =0.05/5;
Para0040.wCf      =0.02;
Para0040.wLc      =0.01;
Para0040.Rc       =0.01/5;
Para0040.Xov      =0.01;
Para0040.N        =1;       % (Hz) matching ratio
Para0040.R        =0.05;    % (Hz) virtual damping resistance
Para0040.fvdq     =300;   % (Hz) vdq bandwidth
Para0040.fidq     =600;  % (Hz) current control bandwidth
Para0040.w0       =W0;
Para0040.Sair     =1000;  % (𝑊/𝑚2) PV illumination intensity
Para0040.Tair     =25;    % (°C) PV temperature
Para0040.C_dc     =1;
Para0040.v_pv_ref =0.8;   % pv output voltage(MPPT) 
Para0040.v_dc_ref =1;   % dc voltage 
Para0040.fvdc     =100; % (Hz) vdc bandwidth
Para0040.fidc     =500; % (Hz) ibat bandwidth
 

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
Para2000.R         = 0.05;
Para2000.K         = 1;
Para2000.N         = 1;
Para2000.v_dc_ref  = 1;

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
            ParaCell{i} = Para0030;     % Battery Energy Storage System
        case 4
            ParaCell{i} = Para0040;     % Photovoltaic
            % PV
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

% Replace the default data by customized data
%
% Notes: 
% This method can reduce the calculation time of "for" loop. The "for" loop
% runs only when "row" is not empty.
%
% The sequence of cases are determined by the excel form. This also
% decouples the sequence between the excel form and the system object.
for i = 1:length(row)
  	AppBus   = AppBusCell{row(i)};
	AppType	= AppTypeCell{row(i)};
 	UserValue 	= ListApp(row(i),column(i));     % Customized value
    SwitchFlag = column(i)-2;                   	% Find the updated parameter
  	if floor(AppType/10) == 0                    % Synchronous machine
        switch SwitchFlag 
         	case 1; ParaCell{row(i)}.J  = UserValue;
            case 2; ParaCell{row(i)}.D  = UserValue;
            case 3; ParaCell{row(i)}.wL = UserValue;
            case 4; ParaCell{row(i)}.R  = UserValue; 
            otherwise
                error(['Error: paramter overflow, bus ' num2str(AppBus) 'type ' num2str(AppType) '.']);
        end
    elseif (floor(AppType/10) == 1)              % Grid-following inverter
        switch SwitchFlag
            case 1; ParaCell{row(i)}.V_dc   = UserValue;
            case 2; ParaCell{row(i)}.C_dc   = UserValue;
            case 3; ParaCell{row(i)}.wLf    = UserValue;
            case 4; ParaCell{row(i)}.R      = UserValue;
            case 5; ParaCell{row(i)}.f_v_dc = UserValue;
            case 6; ParaCell{row(i)}.f_pll  = UserValue;
            case 7; ParaCell{row(i)}.f_i_dq = UserValue;
            otherwise
                error(['Error: parameter overflow, bus ' num2str(AppBus) 'type ' num2str(AppType) '.']);
        end
    elseif floor(AppType/10) == 2                % Grid-forming inverter
        switch SwitchFlag
            case 1;  ParaCell{row(i)}.wLf     = UserValue;
          	case 2;  ParaCell{row(i)}.Rf      = UserValue;
          	case 3;  ParaCell{row(i)}.wCf     = UserValue;
           	case 4;  ParaCell{row(i)}.wLc  	  = UserValue;
         	case 5;  ParaCell{row(i)}.Rc  	  = UserValue;
           	case 6;  ParaCell{row(i)}.Xov 	  = UserValue;
            case 7;  ParaCell{row(i)}.Dw      = UserValue;
            case 8;  ParaCell{row(i)}.fdroop  = UserValue;
          	case 9;  ParaCell{row(i)}.fvdq    = UserValue;
          	case 10; ParaCell{row(i)}.fidq    = UserValue; 
            otherwise
                error(['Error: parameter overflow, bus ' num2str(AppBus) 'type ' num2str(AppType) '.']);
        end
    elseif floor(AppType/10) == 3                % Battery
        switch SwitchFlag
            case 1;  ParaCell{row(i)}.wLf      = UserValue;
            case 2;  ParaCell{row(i)}.Rf       = UserValue;
            case 3;  ParaCell{row(i)}.wCf      = UserValue;
            case 4;  ParaCell{row(i)}.wLc  	   = UserValue;
            case 5;  ParaCell{row(i)}.Rc  	   = UserValue;
            case 6;  ParaCell{row(i)}.Xov 	   = UserValue;
            case 7;  ParaCell{row(i)}.Dw       = UserValue;
            case 8;  ParaCell{row(i)}.fdroop   = UserValue;
            case 9;  ParaCell{row(i)}.fvdq     = UserValue;
            case 10; ParaCell{row(i)}.fidq     = UserValue; 
            case 11; ParaCell{row(i)}.L_dc     = UserValue;
            case 12; ParaCell{row(i)}.C_dc     = UserValue;   
            case 13; ParaCell{row(i)}.v_ocv    = UserValue;
            case 14; ParaCell{row(i)}.R_bat    = UserValue;    
            case 15; ParaCell{row(i)}.v_dc_ref = UserValue;    
            case 16; ParaCell{row(i)}.fvdc     = UserValue; 
            case 17; ParaCell{row(i)}.fibat    = UserValue;  
            otherwise
                error(['Error: parameter overflow, bus ' num2str(AppBus) 'type ' num2str(AppType) '.']);
        end
    elseif floor(AppType/10) == 4                % PV
        switch SwitchFlag
            case 1;  ParaCell{row(i)}.wLf      = UserValue;
            case 2;  ParaCell{row(i)}.Rf       = UserValue;
            case 3;  ParaCell{row(i)}.wCf      = UserValue;
            case 4;  ParaCell{row(i)}.wLc  	   = UserValue;
            case 5;  ParaCell{row(i)}.Rc  	   = UserValue;
            case 6;  ParaCell{row(i)}.Xov 	   = UserValue;
            case 7;  ParaCell{row(i)}.N        = UserValue;
            case 8;  ParaCell{row(i)}.R        = UserValue;
            case 9;  ParaCell{row(i)}.fvdq     = UserValue;
            case 10; ParaCell{row(i)}.fidq     = UserValue; 
            case 11; ParaCell{row(i)}.Sair     = UserValue;
            case 12; ParaCell{row(i)}.Tair     = UserValue;   
            case 13; ParaCell{row(i)}.C_dc     = UserValue;
            case 14; ParaCell{row(i)}.v_pv_ref = UserValue;    
            case 15; ParaCell{row(i)}.v_dc_ref = UserValue;    
            case 16; ParaCell{row(i)}.fvdc     = UserValue; 
            case 17; ParaCell{row(i)}.fidc     = UserValue;  
            otherwise
                error(['Error: parameter overflow, bus ' num2str(AppBus) 'type ' num2str(AppType) '.']);
        end
    elseif floor(AppType/10) == 101 % Grid-feeding buck
        switch SwitchFlag
            case 1;  ParaCell{row(i)}.Vdc   = UserValue;
          	case 2;  ParaCell{row(i)}.Cdc   = UserValue;
          	case 3;  ParaCell{row(i)}.wL    = UserValue;
           	case 4;  ParaCell{row(i)}.R  	= UserValue;
         	case 5;  ParaCell{row(i)}.fi  	= UserValue;
           	case 6;  ParaCell{row(i)}.fvdc 	= UserValue;
            otherwise
                error(['Error: parameter overflow, bus ' num2str(AppBus) 'type ' num2str(AppType) '.']);
        end
    elseif floor(AppType/10) == 200 % Interlink ac-dc converter
        switch SwitchFlag
            case 1;  ParaCell{row(i)}.C_dc     = UserValue;
            case 2;  ParaCell{row(i)}.wL_ac    = UserValue;
            case 3;  ParaCell{row(i)}.R_ac     = UserValue;
            case 4;  ParaCell{row(i)}.wL_dc    = UserValue;
            case 5;  ParaCell{row(i)}.R_dc     = UserValue;
            case 6;  ParaCell{row(i)}.fidq     = UserValue;
            case 7;  ParaCell{row(i)}.fvdc     = UserValue;
            case 8;  ParaCell{row(i)}.fpll     = UserValue;
            case 9;  ParaCell{row(i)}.R        = UserValue; 
            case 10; ParaCell{row(i)}.K        = UserValue;
            case 11; ParaCell{row(i)}.N        = UserValue;
            case 12; ParaCell{row(i)}.v_dc_ref = UserValue;
            otherwise
                error(['Error: parameter overflow, bus ' num2str(AppBus) 'type ' num2str(AppType) '.']);
        end
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