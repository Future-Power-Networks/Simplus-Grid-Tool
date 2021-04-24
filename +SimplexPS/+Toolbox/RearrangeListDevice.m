% This function re-arranges the netlist data of devices.

% Author(s): Yitong Li, Yunjie Gu

%% Notes
%
% The device model is in load convention.

function [DeviceBusCell,DeviceTypeCell,ParaCell,N_Device] = RearrangeListDevice(UserData,W0,ListBus)

%% Load data
[ListDevice,ListDeviceChar]	 = xlsread(UserData,'Device');

%% Rearrange data
[N_Device,ColumnMax_Device] = size(ListDevice);
ListDeviceBus = ListDevice(:,1);
ListDeviceType = ListDevice(:,2);
ListDeviceBusChar = ListDeviceChar(:,1);

for k1 = 1:length(ListDeviceBusChar)
    if strcmpi(ListDeviceBusChar{k1},'Bus No.')
        break;
    end
end
ListDeviceBusChar = ListDeviceBusChar(k1+1:end);

% Get the device bus in cell form
for n = 1:N_Device
    if ~isnan(ListDeviceBus(n))
        DeviceBusCell{n} = ListDeviceBus(n);
    else
        DeviceBusCell{n} = str2num(ListDeviceBusChar{n});
        [~,~,AreaType]= SimplexPS.Toolbox.CheckBus(DeviceBusCell{n}(1),ListBus);
        if AreaType == 2 % If the first bus is dc bus, then swap
            [DeviceBusCell{n}(1),DeviceBusCell{n}(2)] = deal(DeviceBusCell{n}(2),DeviceBusCell{n}(1));
        end
    end
end

% Get the device type in cell form
for n = 1:N_Device
    DeviceTypeCell{n} = ListDevice(n,2);
end

% Re-order the device sequence
% ListDevice = sortrows(ListDevice,1);
% No re-order for now.

% Error check
if (ColumnMax_Device>12)
    error(['Error: Device data overflow.']); 
end

[~,ModeBus] = SimplexPS.CellMode(DeviceBusCell);
if ModeBus~=1
    error(['Error: For each bus, one and only one device has to be connected.']);
end

%% Default AC device data
% ======================================
% Synchronous generator
% ======================================
Para0000.J  = 3.5;
Para0000.D  = 1;
Para0000.L  = 0.1;
Para0000.R  = 0.01;
Para0000.w0 = W0;

% ======================================
% Grid-following VSI (PLL-controlled)
% ======================================
% Dc link
Para0010.V_dc   = 2.5;
Para0010.C_dc   = 2*0.1*Para0010.V_dc^2;
Para0010.f_v_dc = 20;           % (Hz) bandwidth, vdc

% Ac filter
Para0010.wLf = 0.03;
Para0010.R   = 0.01;

% PLL
Para0010.f_pll      = 20;       % (Hz) bandwidth, PLL
Para0010.f_tau_pll  = 200;      % (Hz) bandwidth, PLL low pass filter

% Current loop
Para0010.f_i_dq = 500;      	% (Hz) bandwidth, idq
Para0010.w0 = W0;   

% ======================================
% Grid-forming VSI (Droop-Controlled)
% ======================================
% Bandwidth
Para0020.w_i_ldq  = 500*2*pi;     % (rad/s), current loop bandwidth
Para0020.w_v_odq  = 250*2*pi;     % (rad/s), voltage loop bandwidth 
Para0020.wf       = 20*2*pi;      % (rad/s), droop filter bandwidth

% Ac filter
Para0020.Lf = 0.05/W0;
Para0020.Rf = 0.05/5;
Para0020.Cf = 0.02/W0;
Para0020.Lc = 0.01/W0;
Para0020.Rc = 0.01/5;

% Droop loop
Para0020.Dw = 5/100*W0;
% Para20.Dv = 5/100;
% Para20.P0 = 0;
% Para20.Q0 = 0;
Para0020.w0 = W0;

% Voltage and current loop
% Para20.Gi_cd = 0;
% Para20.Gv_cd = 0;
% Para20.Fv = 0;
% Para20.Fi = 0;

% Virtual impedance
% Para20.Rov = 0;
Para0020.Xov = 0;

% ======================================
% Ac infinite bus (short-circuit in small-signal)
% ======================================
Para0090 = [];

% ======================================
% Ac floating bus (open-circuit)
% ======================================
Para0100 = [];

%% Default DC device data

% ======================================
% Grid-feeding buck
% ======================================
% Bandwidth
w_vdc	= 10*2*pi; 	% (rad/s) bandwidth, vdc
w_i     = 500*2*pi;	% (rad/s) bandwidth, i

% Dc link
Para1010.V_dc   	= 2;
Para1010.C_dc       = 2*0.1*Para1010.V_dc^2;
Para1010.kp_v_dc	= Para1010.V_dc*Para1010.C_dc*w_vdc;
Para1010.ki_v_dc	= Para1010.kp_v_dc*w_vdc/4;

% Dc-grid-side filter
Para1010.L        = 0.05/W0;
Para1010.R        = 0.01;

% Current loop
Para1010.kp_i  = Para1010.L * w_i;         % P
Para1010.ki_i  = Para1010.L * w_i^2 /4;    % I

% ======================================
% Dc infinite bus (short-circuit in small-signal)
% ======================================
Para1090 = [];

% ======================================
% Dc floating bus (open-circuit)
% ======================================
Para1100 = [];

%% Default hybrid device data

% ======================================
% Interlink ac-dc converter
% ======================================

% Bandwidth
w_vdc     = 10*2*pi; 	% (rad/s) bandwidth, vdc
w_pll     = 10*2*pi;  	% (rad/s) bandwidth, pll
w_i     = 500*2*pi; 	% (rad/s) bandwidth, idq
w_tau_pll = 200*2*pi;   % (rad/s) PLL filter bandwidth

% DC link loop
V_dc = 1;
Para2000.kp_v_dc = V_dc*Para0010.C_dc*w_vdc;
Para2000.ki_v_dc = Para2000.kp_v_dc*w_vdc/4;

% Dc filter
Para2000.C_dc 	= 2*0.1*V_dc^2 * 8;
Para2000.L_dc 	= 0.01/W0;
Para2000.R_dc	= 0.01/5;

% Ac filter
Para2000.L_ac 	= 0.05/W0;
Para2000.R_ac 	= 0.01;

% PLL
Para2000.kp_pll   = w_pll;
Para2000.ki_pll   = Para2000.kp_pll * w_pll/4; 
Para2000.tau_pll  = 1/w_tau_pll;

% Current loop
Para2000.kp_i_dq  = Para2000.L_ac * w_i;         % P
Para2000.ki_i_dq  = Para2000.L_ac * w_i^2 /4;    % I
Para2000.w0       = W0;   

%% Re-arrange device data
% Get the size of netlist
[N_Device,ColumnMax_Device] = size(ListDevice);

% Find the index of user-defined data
netlist_device_NaN = isnan(ListDevice(:,3:ColumnMax_Device));
[row,column] = find(netlist_device_NaN == 0);     
column = column+2;

% Initialize the device parameters by default parameters
for i = 1:N_Device
    DeviceBus   = DeviceBusCell{i};
    DeviceType  = ListDeviceType(i);
    switch floor(DeviceType/10)
        % ### AC devices
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
            ParaCell{i} = Para0100;     % Ac floating bus, i.e., no device
        
        % ### DC devices
        case 101
            ParaCell{i} = Para1010;     % Grid-following buck
        case 109
            ParaCell{i} = Para1090;     % Dc infinite bus
        case 110
            ParaCell{i} = Para1100;     % Ac floating bus, i.e., no device
            
        % ### Hybrid ac-dc devices
        case 200
            ParaCell{i} = Para2000;     % Interlinking ac-dc converter
            
        % ### Error check
        otherwise
            error(['Error: device type, bus ' num2str(DeviceBus) ' type ' num2str(DeviceType) '.']);
    end
end

% Replace the default data by customized data
% Notes: 
% This method can reduce the calculation time of "for" loop.
% The "for" loop runs only when "row" is not empty.
%
% The sequence of cases are determined by the excel form. This also
% decouples the sequence between the excel form and the system object.
for i = 1:length(row)
  	DeviceBus   = DeviceBusCell{row(i)};
	DeviceType	= ListDeviceType(row(i));
 	UserValue 	= ListDevice(row(i),column(i));      % Customized value
    SwitchFlag = column(i)-2;                         	% Find the updated parameter
  	if floor(DeviceType/10) == 0
        switch SwitchFlag 
         	case 1; ParaCell{row(i)}.J  = UserValue;
            case 2; ParaCell{row(i)}.D  = UserValue;
            case 3; ParaCell{row(i)}.wL = UserValue;
            case 4; ParaCell{row(i)}.R  = UserValue; 
            otherwise
                error(['Error: paramter overflow, bus ' num2str(DeviceBus) 'type ' num2str(DeviceType) '.']);
        end
    elseif (floor(DeviceType/10) == 1)
        switch SwitchFlag
            case 1; ParaCell{row(i)}.V_dc   = UserValue;
            case 2; ParaCell{row(i)}.C_dc   = UserValue;
            case 3; ParaCell{row(i)}.wL     = UserValue;
            case 4; ParaCell{row(i)}.R      = UserValue;
            case 5; ParaCell{row(i)}.f_v_dc = UserValue;
            case 6; ParaCell{row(i)}.f_pll  = UserValue;
            case 7; ParaCell{row(i)}.f_i_dq = UserValue;
            otherwise
                error(['Error: parameter overflow, bus ' num2str(DeviceBus) 'type ' num2str(DeviceType) '.']);
        end
    elseif floor(DeviceType/10) == 2
        switch SwitchFlag
            case 1;  ParaCell{row(i)}.Lf      = UserValue/W0;
          	case 2;  ParaCell{row(i)}.Rf      = UserValue;
          	case 3;  ParaCell{row(i)}.Cf      = UserValue/W0;
           	case 4;  ParaCell{row(i)}.Lc  	  = UserValue/W0;
         	case 5;  ParaCell{row(i)}.Rc  	  = UserValue;
           	case 6;  ParaCell{row(i)}.Xov 	  = UserValue;
            case 7;  ParaCell{row(i)}.Dw      = UserValue*W0;
            case 8;  ParaCell{row(i)}.wf      = UserValue*2*pi;
          	case 9;  ParaCell{row(i)}.w_v_odq = UserValue*2*pi;
          	case 10; ParaCell{row(i)}.w_i_ldq = UserValue*2*pi;
            otherwise
                error(['Error: parameter overflow, bus ' num2str(DeviceBus) 'type ' num2str(DeviceType) '.']);
        end
    end
end

end