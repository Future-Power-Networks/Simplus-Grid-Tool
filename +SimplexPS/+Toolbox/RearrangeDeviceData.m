% This function re-arranges the device data

% Author(s): Yitong Li, Yunjie Gu

%% Notes
%
% The device model is in load convention.

function [CellDeviceType,CellPara] = RearrangeDeviceData(NetlistDevice,W0)

%% Default device data
% ======================================
% Synchronous generator
% ======================================
Para00.J  = 3.5*2/W0^2;
Para00.D  = 1/W0^2;
Para00.L  = 0.1/W0;
Para00.R  = 0.01;
Para00.w0 = W0;

% ======================================
% Grid-following VSI (PLL-controlled)
% ======================================
% Bandwidth
w_vdc     = 20*2*pi; 	% (rad/s) bandwidth, vdc
w_pll     = 20*2*pi;  	% (rad/s) bandwidth, pll
w_idq     = 500*2*pi; 	% (rad/s) bandwidth, idq
w_tau_pll = 200*2*pi;   % (rad/s) PLL filter bandwidth

% DC link
Para10.V_dc   	= 2.5;
Para10.C_dc 	= 2*0.1*Para10.V_dc^2;
Para10.kp_v_dc	= Para10.V_dc*Para10.C_dc*w_vdc;
Para10.ki_v_dc	= Para10.kp_v_dc*w_vdc/4;

% AC filter
Para10.L        = 0.03/W0;
Para10.R        = 0.01;

% PLL
Para10.kp_pll   = w_pll;
Para10.ki_pll   = Para10.kp_pll * w_pll/4; 
Para10.tau_pll  = 1/w_tau_pll;

% Current loop
Para10.k_pf     = 0;
Para10.kp_i_dq  = Para10.L * w_idq;         % P
Para10.ki_i_dq  = Para10.L * w_idq^2 /4;    % I
Para10.w0       = W0;   
Para10.Gi_cd    = 0;                        % Cross-decoupling gain

% ======================================
% Grid-forming VSI (Droop-Controlled)
% ======================================
% Bandwidth
Para20.w_i_ldq  = 500*2*pi;     % (rad/s), current loop bandwidth
Para20.w_v_odq  = 250*2*pi;     % (rad/s), voltage loop bandwidth 
Para20.wf       = 20*2*pi;      % (rad/s), droop filter bandwidth

% DC link
% Para10.V_dc	= 2.5;
% Para10.C_dc	= 2*0.1*Para10.V_dc^2;

% AC filter
Para20.Lf = 0.05/W0;
Para20.Rf = 0.05/5;
Para20.Cf = 0.02/W0;
Para20.Lc = 0.01/W0;
Para20.Rc = 0.01/5;

% Droop loop
Para20.Dw = 5/100*W0;
% Para20.Dv = 5/100;
% Para20.P0 = 0;
% Para20.Q0 = 0;
Para20.w0 = W0;

% Voltage and current loop
% Para20.Gi_cd = 0;
% Para20.Gv_cd = 0;
% Para20.Fv = 0;
% Para20.Fi = 0;

% Virtual impedance
% Para20.Rov = 0;
Para20.Xov = 0;

% ======================================
% Infinite bus (short-circuit in small-signal)
% ======================================
Para90 = [];

% ======================================
% Floating bus (open-circuit)
% ======================================
Para100 = [];

%% Re-arrange device data
% Get the size of netlist
[N_Device,ColumnMax_Device] = size(NetlistDevice);

% Find the index of user-defined data
netlist_device_NaN = isnan(NetlistDevice(:,3:ColumnMax_Device));
[row,column] = find(netlist_device_NaN == 0);     
column = column+2;

% Initialize the device parameters by default parameters
for i = 1:N_Device
    bus   = NetlistDevice(i,1);
    type  = NetlistDevice(i,2);
    CellDeviceType{i} = type;
    switch floor(type/10)
        case 0     
            CellPara{i} = Para00;   % Synchronous machine
        case 1
            CellPara{i} = Para10;   % Grid-following inverter
         case 2
            CellPara{i} = Para20;   % Grid-forming inverter
        case 9
            CellPara{i} = Para90;   % Inifnite bus
        case 10
            CellPara{i} = Para100;  % Floating bus, i.e., no device
        otherwise
            error(['Error: device type, bus ' num2str(bus) ' type ' num2str(type) '.']);
    end
end

% Replace the default data by customized data
% Notes: 
% This method can reduce the calculation time of "for" loop.
% The "for" loop runs only when "row" is not empty.
for i = 1:length(row)
  	bus         = NetlistDevice(row(i),1);
	type        = NetlistDevice(row(i),2);
 	user_value 	= NetlistDevice(row(i),column(i));      % Customized value
    switch_flag = column(i)-2;                         	% Find the updated parameter
  	if floor(type/10) == 0
        switch switch_flag 
         	case 1; CellPara{row(i)}.J = user_value*2/W0^2;
            case 2; CellPara{row(i)}.D = user_value/W0^2;
            case 3; CellPara{row(i)}.L = user_value/W0;
            case 4; CellPara{row(i)}.R = user_value; 
            otherwise
                error(['Error: paramter overflow, bus ' num2str(bus) 'type ' num2str(type) '.']);
        end
    elseif floor(type/10) == 1
        switch switch_flag
            case 1; CellPara{row(i)}.V_dc     = user_value;
            case 2; CellPara{row(i)}.C_dc     = user_value;
            case 3; CellPara{row(i)}.L        = user_value/W0;
            case 4; CellPara{row(i)}.R        = user_value;
            case 5; CellPara{row(i)}.kp_v_dc  = CellPara{row(i)}.V_dc*CellPara{row(i)}.C_dc*(user_value*2*pi);
                    CellPara{row(i)}.ki_v_dc  = CellPara{row(i)}.kp_v_dc*(user_value*2*pi)/4;
            case 6; CellPara{row(i)}.kp_pll   = user_value*2*pi;
                    CellPara{row(i)}.ki_pll   = CellPara{row(i)}.kp_pll*(user_value*2*pi)/4; 
            case 7; CellPara{row(i)}.kp_i_dq  = CellPara{row(i)}.L*(user_value*2*pi);
                    CellPara{row(i)}.ki_i_dq  = CellPara{row(i)}.kp_i_dq *(user_value*2*pi)/4;
            otherwise
                error(['Error: parameter overflow, bus ' num2str(bus) 'type ' num2str(type) '.']);
        end
    elseif floor(type/10) == 2
        switch switch_flag
            case 1;  CellPara{row(i)}.Lf      = user_value/W0;
          	case 2;  CellPara{row(i)}.Rf      = user_value;
          	case 3;  CellPara{row(i)}.Cf      = user_value/W0;
           	case 4;  CellPara{row(i)}.Lc  	  = user_value/W0;
         	case 5;  CellPara{row(i)}.Rc  	  = user_value/W0;
           	case 6;  CellPara{row(i)}.Xov 	  = user_value;
            case 7;  CellPara{row(i)}.Dw      = user_value*W0;
            case 8;  CellPara{row(i)}.wf      = user_value*2*pi;
          	case 9;  CellPara{row(i)}.w_v_odq = user_value*2*pi;
          	case 10; CellPara{row(i)}.w_i_ldq = user_value*2*pi;
            otherwise
                error(['Error: parameter overflow, bus ' num2str(bus) 'type ' num2str(type) '.']);
        end
    end
end

end