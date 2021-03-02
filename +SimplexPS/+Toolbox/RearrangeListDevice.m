% This function re-arranges the netlist data of devices.

% Author(s): Yitong Li, Yunjie Gu

%% Notes
%
% The device model is in load convention.

function [CellDeviceType,CellPara,N_Device] = RearrangeListDevice(Name_Netlist,W0)

%% Load data
ListDevice	 = xlsread(Name_Netlist,'Device');

%% Rearrange data
% Re-order the device sequence
ListDevice = sortrows(ListDevice,1);

% Error check
[N_Device,ColumnMax_Device] = size(ListDevice);
if (ColumnMax_Device>12)
    error(['Error: Device data overflow.']); 
end
[~,ModeBus] = mode(ListDevice(:,1));
if ModeBus~=1
    error(['Error: For each bus, one and only one device has to be connected.']);
end

%% Default AC device data
% ======================================
% Synchronous generator
% ======================================
Para0000.J  = 3.5*2/W0^2;
Para0000.D  = 1/W0^2;
Para0000.L  = 0.1/W0;
Para0000.R  = 0.01;
Para0000.w0 = W0;

% ======================================
% Grid-following VSI (PLL-controlled)
% ======================================
% Bandwidth
w_vdc     = 20*2*pi; 	% (rad/s) bandwidth, vdc
w_pll     = 20*2*pi;  	% (rad/s) bandwidth, pll
w_idq     = 500*2*pi; 	% (rad/s) bandwidth, idq
w_tau_pll = 200*2*pi;   % (rad/s) PLL filter bandwidth

% DC link
Para0010.V_dc   	= 2.5;
Para0010.C_dc       = 2*0.1*Para0010.V_dc^2;
Para0010.kp_v_dc	= Para0010.V_dc*Para0010.C_dc*w_vdc;
Para0010.ki_v_dc	= Para0010.kp_v_dc*w_vdc/4;

% AC filter
Para0010.L        = 0.03/W0;
Para0010.R        = 0.01;

% PLL
Para0010.kp_pll   = w_pll;
Para0010.ki_pll   = Para0010.kp_pll * w_pll/4; 
Para0010.tau_pll  = 1/w_tau_pll;

% Current loop
Para0010.k_pf     = 0;
Para0010.kp_i_dq  = Para0010.L * w_idq;         % P
Para0010.ki_i_dq  = Para0010.L * w_idq^2 /4;    % I
Para0010.w0       = W0;   
Para0010.Gi_cd    = 0;                        % Cross-decoupling gain

% ======================================
% Grid-forming VSI (Droop-Controlled)
% ======================================
% Bandwidth
Para0020.w_i_ldq  = 500*2*pi;     % (rad/s), current loop bandwidth
Para0020.w_v_odq  = 250*2*pi;     % (rad/s), voltage loop bandwidth 
Para0020.wf       = 20*2*pi;      % (rad/s), droop filter bandwidth

% DC link
% Para10.V_dc	= 2.5;
% Para10.C_dc	= 2*0.1*Para10.V_dc^2;

% AC filter
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
% AC infinite bus (short-circuit in small-signal)
% ======================================
Para0090 = [];

% ======================================
% AC floating bus (open-circuit)
% ======================================
Para0100 = [];

%% Default DC device data

% ======================================
% Grid-following buck
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
Para1010.kp_i  = Para0010.L * w_i;         % P
Para1010.ki_i  = Para0010.L * w_i^2 /4;    % I

% ======================================
% DC infinite bus (short-circuit in small-signal)
% ======================================
Para1090 = [];

% ======================================
% DC floating bus (open-circuit)
% ======================================
Para1100 = [];

%% Default hybrid device data

%% Re-arrange device data
% Get the size of netlist
[N_Device,ColumnMax_Device] = size(ListDevice);

% Find the index of user-defined data
netlist_device_NaN = isnan(ListDevice(:,3:ColumnMax_Device));
[row,column] = find(netlist_device_NaN == 0);     
column = column+2;

% Initialize the device parameters by default parameters
for i = 1:N_Device
    bus   = ListDevice(i,1);
    type  = ListDevice(i,2);
    CellDeviceType{i} = type;
    switch floor(type/10)
        case 0     
            CellPara{i} = Para0000;   % Synchronous machine
        case 1
            CellPara{i} = Para0010;   % Grid-following inverter
         case 2
            CellPara{i} = Para0020;   % Grid-forming inverter
        case 9
            CellPara{i} = Para0090;   % Inifnite bus
        case 10
            CellPara{i} = Para0100;  % Floating bus, i.e., no device
        otherwise
            error(['Error: device type, bus ' num2str(bus) ' type ' num2str(type) '.']);
    end
end

% Replace the default data by customized data
% Notes: 
% This method can reduce the calculation time of "for" loop.
% The "for" loop runs only when "row" is not empty.
for i = 1:length(row)
  	bus         = ListDevice(row(i),1);
	type        = ListDevice(row(i),2);
 	UserValue 	= ListDevice(row(i),column(i));      % Customized value
    SwitchFlag = column(i)-2;                         	% Find the updated parameter
  	if floor(type/10) == 0
        switch SwitchFlag 
         	case 1; CellPara{row(i)}.J = UserValue*2/W0^2;
            case 2; CellPara{row(i)}.D = UserValue/W0^2;
            case 3; CellPara{row(i)}.L = UserValue/W0;
            case 4; CellPara{row(i)}.R = UserValue; 
            otherwise
                error(['Error: paramter overflow, bus ' num2str(bus) 'type ' num2str(type) '.']);
        end
    elseif floor(type/10) == 1
        switch SwitchFlag
            case 1; CellPara{row(i)}.V_dc     = UserValue;
            case 2; CellPara{row(i)}.C_dc     = UserValue;
            case 3; CellPara{row(i)}.L        = UserValue/W0;
            case 4; CellPara{row(i)}.R        = UserValue;
            case 5; CellPara{row(i)}.kp_v_dc  = CellPara{row(i)}.V_dc*CellPara{row(i)}.C_dc*(UserValue*2*pi);
                    CellPara{row(i)}.ki_v_dc  = CellPara{row(i)}.kp_v_dc*(UserValue*2*pi)/4;
            case 6; CellPara{row(i)}.kp_pll   = UserValue*2*pi;
                    CellPara{row(i)}.ki_pll   = CellPara{row(i)}.kp_pll*(UserValue*2*pi)/4; 
            case 7; CellPara{row(i)}.kp_i_dq  = CellPara{row(i)}.L*(UserValue*2*pi);
                    CellPara{row(i)}.ki_i_dq  = CellPara{row(i)}.kp_i_dq *(UserValue*2*pi)/4;
            otherwise
                error(['Error: parameter overflow, bus ' num2str(bus) 'type ' num2str(type) '.']);
        end
    elseif floor(type/10) == 2
        switch SwitchFlag
            case 1;  CellPara{row(i)}.Lf      = UserValue/W0;
          	case 2;  CellPara{row(i)}.Rf      = UserValue;
          	case 3;  CellPara{row(i)}.Cf      = UserValue/W0;
           	case 4;  CellPara{row(i)}.Lc  	  = UserValue/W0;
         	case 5;  CellPara{row(i)}.Rc  	  = UserValue/W0;
           	case 6;  CellPara{row(i)}.Xov 	  = UserValue;
            case 7;  CellPara{row(i)}.Dw      = UserValue*W0;
            case 8;  CellPara{row(i)}.wf      = UserValue*2*pi;
          	case 9;  CellPara{row(i)}.w_v_odq = UserValue*2*pi;
          	case 10; CellPara{row(i)}.w_i_ldq = UserValue*2*pi;
            otherwise
                error(['Error: parameter overflow, bus ' num2str(bus) 'type ' num2str(type) '.']);
        end
    end
end

end