% This function re-arranges the device data

% Author(s): Yitong Li, Yunjie Gu, Yue Zhu

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
% Synchronous generato ----- Full Model
% ======================================
Para30.X=0.0125;
Para30.R=0;
Para30.Xd=0.1; %synchronous reactance in d axis
Para30.Xd1=0.031; %transient reactance
Para30.Xd2=0.025; %subtransient reactance
Para30.Td1=10.2; %d-axis open circuit transient time constant
Para30.Td2=0.05; %d-axis open circuit sub-transient time constant
Para30.Xq=0.069;
Para30.Xq1=0.0416667;
Para30.Xq2=0.025;
Para30.Tq1=1.5;
Para30.Tq2=0.035;
Para30.H=42;
Para30.D=0;
Para30.TR=0.01;
Para30.KA=1;
Para30.TA=0.02;
Para30.VRmax=10;
Para30.VRmin=-10;
Para30.KE=1;
Para30.TE=0.785;
Para30.E1=3.9267;
Para30.SE1=0.07;
Para30.E2=5.2356;
Para30.SE2=0.91;
Para30.KF=0.03;
Para30.TF=1;
Para30.KP=200;
Para30.KI=50;
Para30.KD=50;
Para30.TD=0.01;
Para30.KPSS=20;
Para30.TW=15;
Para30.T11=0.15;
Para30.T12=0.04;
Para30.T21=0.15;
Para30.T22=0.04;
Para30.T31=0.15;
Para30.T32=0.04;
Para30.VSSmax=0.2;
Para30.VSSmin=-0.05;
Para30.Rgov=0.05;
Para30.T1gov=0.8;
Para30.T2gov=1.8;
Para30.T3gov=6;
Para30.Dtgov=0.2;

% ======================================
% Grid-following VSI (PLL-controlled)
% ======================================
% Bandwidth
% w_vdc     = 20*2*pi; 	% (rad/s) bandwidth, vdc
% w_pll     = 20*2*pi;  	% (rad/s) bandwidth, pll
% w_idq     = 500*2*pi; 	% (rad/s) bandwidth, idq
% w_tau_pll = 200*2*pi;   % (rad/s) PLL filter bandwidth

% DC link
% Para10.V_dc   	= 2.5;
% Para10.C_dc 	= 2*0.1*Para10.V_dc^2;
% Para10.kp_v_dc	= Para10.V_dc*Para10.C_dc*w_vdc;
% Para10.ki_v_dc	= Para10.kp_v_dc*w_vdc/4;
% 
% % AC filter
% Para10.L        = 0.03/W0;
% Para10.R        = 0.01;
% 
% % PLL
% Para10.kp_pll   = w_pll;
% Para10.ki_pll   = Para10.kp_pll * w_pll/4; 
% Para10.tau_pll  = 1/w_tau_pll;
% 
% % Current loop
% Para10.k_pf     = 0;
% Para10.kp_i_dq  = Para10.L * w_idq;         % P
% Para10.ki_i_dq  = Para10.L * w_idq^2 /4;    % I
% Para10.w0       = W0;   
% Para10.Gi_cd    = 0;                        % Cross-decoupling gain
Para10.Vdc=2.5;
Para10.Cdc=1.25;
Para10.wL=0.03;
Para10.R = 0.01;
Para10.fvdc=20;
Para10.fpll=20;
Para10.fidq=500;

% ======================================
% Grid-forming VSI (Droop-Controlled)
% ======================================
Para20.wLf=0.05;
Para20.Rf=0.01;
Para20.wCf=0.02;
Para20.wLc=0.01;
Para20.Rc=0.002;
Para20.Xov=0;
Para20.Dw=0.002;
Para20.fdroop=20;
Para20.fvdc=250;
Para20.fidq=500;

%                             Para20.w_i_ldq  = 500*2*pi;     % (rad/s), current loop bandwidth
%                             Para20.w_v_odq  = 250*2*pi;     % (rad/s), voltage loop bandwidth 
%                             Para20.wf       = 20*2*pi;      % (rad/s), droop filter bandwidth
%                             Para20.Lf = 0.05/W0;
%                             Para20.Rf = 0.05/5;
%                             Para20.Cf = 0.02/W0;
%                             Para20.Lc = 0.01/W0;
%                             Para20.Rc = 0.01/5;
%                             Para20.Dw = 5/100*W0;
%                             Para20.w0 = W0;
%                             Para20.Xov = 0;

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
        case 3
            CellPara{i} = Para30;   % Synchronous machine full model
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
            case 1; CellPara{row(i)}.Vdc=user_value;
            case 2; CellPara{row(i)}.Cdc=user_value;
            case 3; CellPara{row(i)}.wL=user_value;
            case 4; CellPara{row(i)}.R=user_value;
            case 5; CellPara{row(i)}.fvdc=user_value;
            case 6; CellPara{row(i)}.fpll=user_value;
            case 7; CellPara{row(i)}.fidq=user_value;                  
%             case 1; CellPara{row(i)}.V_dc     = user_value;
%             case 2; CellPara{row(i)}.C_dc     = user_value;
%             case 3; CellPara{row(i)}.L        = user_value/W0;
%             case 4; CellPara{row(i)}.R        = user_value;
%             case 5; CellPara{row(i)}.kp_v_dc  = CellPara{row(i)}.V_dc*CellPara{row(i)}.C_dc*(user_value*2*pi);
%                     CellPara{row(i)}.ki_v_dc  = CellPara{row(i)}.kp_v_dc*(user_value*2*pi)/4;
%             case 6; CellPara{row(i)}.kp_pll   = user_value*2*pi;
%                     CellPara{row(i)}.ki_pll   = CellPara{row(i)}.kp_pll*(user_value*2*pi)/4; 
%             case 7; CellPara{row(i)}.kp_i_dq  = CellPara{row(i)}.L*(user_value*2*pi);
%                     CellPara{row(i)}.ki_i_dq  = CellPara{row(i)}.kp_i_dq *(user_value*2*pi)/4;
            otherwise
                error(['Error: parameter overflow, bus ' num2str(bus) 'type ' num2str(type) '.']);
        end
    elseif floor(type/10) == 2
        switch switch_flag
            case 1;  CellPara{row(i)}.wLf=user_value;
            case 2;  CellPara{row(i)}.Rf=user_value;
            case 3;  CellPara{row(i)}.wCf=user_value;
            case 4;  CellPara{row(i)}.wLc=user_value;
            case 5;  CellPara{row(i)}.Rc=user_value;
            case 6;  CellPara{row(i)}.Xov=user_value;
            case 7;  CellPara{row(i)}.Dw=user_value;
            case 8;  CellPara{row(i)}.fdroop=user_value;
            case 9;  CellPara{row(i)}.fvdc=user_value;
            case 10;  CellPara{row(i)}.fidq=user_value;

%             case 1;  CellPara{row(i)}.Lf      = user_value/W0;
%           	case 2;  CellPara{row(i)}.Rf      = user_value;
%           	case 3;  CellPara{row(i)}.Cf      = user_value/W0;
%            	case 4;  CellPara{row(i)}.Lc  	  = user_value/W0;
%          	case 5;  CellPara{row(i)}.Rc  	  = user_value/W0;
%            	case 6;  CellPara{row(i)}.Xov 	  = user_value;
%             case 7;  CellPara{row(i)}.Dw      = user_value*W0;
%             case 8;  CellPara{row(i)}.wf      = user_value*2*pi;
%           	case 9;  CellPara{row(i)}.w_v_odq = user_value*2*pi;
%           	case 10; CellPara{row(i)}.w_i_ldq = user_value*2*pi;
            otherwise
                error(['Error: parameter overflow, bus ' num2str(bus) 'type ' num2str(type) '.']);
        end
        
    elseif floor(type/10) == 3 %full model
        switch switch_flag
            case 1; CellPara{row(i)}.X=user_value;
            case 2; CellPara{row(i)}.R=user_value;
            case 3; CellPara{row(i)}.Xd=user_value; 
            case 4; CellPara{row(i)}.Xd1=user_value; 
            case 5; CellPara{row(i)}.Xd2=user_value; 
            case 6; CellPara{row(i)}.Td1=user_value; 
            case 7; CellPara{row(i)}.Td2=user_value;
            case 8; CellPara{row(i)}.Xq=user_value;
            case 9; CellPara{row(i)}.Xq1=user_value;
            case 10;CellPara{row(i)}.Xq2=user_value;
            case 11;CellPara{row(i)}.Tq1=user_value;
            case 12;CellPara{row(i)}.Tq2=user_value;
            case 13;CellPara{row(i)}.H=user_value;
            case 14;CellPara{row(i)}.D=user_value;
            case 15;CellPara{row(i)}.TR=user_value;
            case 16;CellPara{row(i)}.KA=user_value;
            case 17;CellPara{row(i)}.TA=user_value;
            case 18;CellPara{row(i)}.VRmax=user_value;
            case 19;CellPara{row(i)}.VRmin=user_value;
            case 20;CellPara{row(i)}.KE=user_value;
            case 21;CellPara{row(i)}.TE=user_value;
            case 22;CellPara{row(i)}.E1=user_value;
            case 23;CellPara{row(i)}.SE1=user_value;
            case 24;CellPara{row(i)}.E2=user_value;
            case 25;CellPara{row(i)}.SE2=user_value;
            case 26;CellPara{row(i)}.KF=user_value;
            case 27;CellPara{row(i)}.TF=user_value;
            case 28;CellPara{row(i)}.KP=user_value;
            case 29;CellPara{row(i)}.KI=user_value;
            case 30;CellPara{row(i)}.KD=user_value;
            case 31;CellPara{row(i)}.TD=user_value;
            case 32;CellPara{row(i)}.KPSS=user_value;
            case 33;CellPara{row(i)}.TW=user_value;
            case 34;CellPara{row(i)}.T11=user_value;
            case 35;CellPara{row(i)}.T12=user_value;
            case 36;CellPara{row(i)}.T21=user_value;
            case 37;CellPara{row(i)}.T22=user_value;
            case 38;CellPara{row(i)}.T31=user_value;
            case 39;CellPara{row(i)}.T32=user_value;
            case 40;CellPara{row(i)}.VSSmax=user_value;
            case 41;CellPara{row(i)}.VSSmin=user_value;
            case 42;CellPara{row(i)}.Rgov=user_value;
            case 43;CellPara{row(i)}.T1gov=user_value;
            case 44;CellPara{row(i)}.T2gov=user_value;
            case 45;CellPara{row(i)}.T3gov=user_value;
            case 46;CellPara{row(i)}.Dtgov=user_value;
            otherwise
                error(['Error: parameter overflow, bus ' num2str(bus) 'type ' num2str(type) '.']);
        end
    end
end

end