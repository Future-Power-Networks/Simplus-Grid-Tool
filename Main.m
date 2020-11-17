% Author(s): Yitong Li, Yunjie Gu

%% 
% Notes:
%
% Please read "README.md" first before using the toolbox.
%
% Please use "Main_Customer.m" rather than this file for running the
% toolbox.

%%
clear all;  % Clear matlab workspace
clc;        % Clear matlab command window
close all;  % Close all figures, etc

fprintf('==================================\n')
fprintf('Start\n')
fprintf('==================================\n')

%%
% ==================================================
% Add folder to path
% ==================================================
% AddFolder2Path();
addpath('GenericFunction');
addpath('Plot');
addpath('SimulinkModel');
addpath('SystemObject');
addpath('Toolbox');
addpath('ForTest');

%%
% ==================================================
% Load customized data
% ==================================================
fprintf('Load customized data.\n')

% ### Load the customized data
% Other available function: readmatrix, csvread ...
Name_Netlist = 'netlist.xlsx';

ListBus    	 = xlsread(Name_Netlist,1);     
ListDevice	 = xlsread(Name_Netlist,2);
ListBasic    = xlsread(Name_Netlist,3);
ListLine     = xlsread(Name_Netlist,4);
ListLineIEEE = xlsread(Name_Netlist,5);
ListAdvance  = xlsread(Name_Netlist,6);

% ### Re-arrange advanced settings
Flag_PowerFlowAlgorithm      	= ListAdvance(5);
Enable_CreateSimulinkModel      = ListAdvance(6);
Enable_PlotPole                 = ListAdvance(7);
Enable_PlotAdmittance           = ListAdvance(8);
Enable_PrintOutput              = ListAdvance(9);
Enable_PlotSwing                = ListAdvance(10);

% ### Re-arrange the simulation data
Fs = ListBasic(1);
Ts = 1/Fs;              % (s), sampling period
Fbase = ListBasic(2);   % (Hz), base frequency
Sbase = ListBasic(3);   % (VA), base power
Vbase = ListBasic(4);   % (V), base voltage
Ibase = Sbase/Vbase;    % (A), base current
Zbase = Sbase/Ibase;    % (Ohm), base impedance
Ybase = 1/Zbase;        % (S), base admittance
Wbase = Fbase*2*pi;     % (rad/s), base angular frequency

% ### Re-arrange the netlist and check error
[ListLine] = RearrangeNetlist_IEEE2Toolbox(ListLine,ListLineIEEE);
[ListBus,ListLine,ListDevice,N_Bus,N_Branch,N_Device] = RearrangeNetlist(ListBus,ListLine,ListDevice);

% ### Re-arrange the device data
[DeviceType,Para] = RearrangeDeviceData(ListDevice,Wbase);

%%
% ==================================================
% Descriptor state space model
% ==================================================

% ### Power flow analysis
fprintf('Do the power flow analysis.\n')
switch Flag_PowerFlowAlgorithm
    case 1  % Gauss-Seidel 
        [PowerFlow,~,~,~,~,~,~,~] = PowerFlow_GS(ListBus,ListLine,Wbase);
    case 2  % Newton-Raphson
        [PowerFlow] = PowerFlow_NR(ListBus,ListLine,Wbase);
    otherwise
        error(['Error: Wrong setting for power flow algorithm.']);
end
ListPowerFlow = RearrangePowerFlow(PowerFlow);
% Move load flow to bus admittance matrix
[ListBus,ListLine,PowerFlow] = Load2SelfBranch(ListBus,ListLine,DeviceType,PowerFlow);
ListPowerFlow_ = RearrangePowerFlow(PowerFlow);

% ### Get the model of lines
fprintf('Get the descriptor state space model of network lines.\n')

[YbusObj,YbusDSS,~] = YbusCalcDSS(ListLine,Wbase);
[~,lsw] = size(YbusDSS.B);
ZbusObj = obj_SwitchInOut(YbusObj,lsw);
[ZbusStateStr,ZbusInputStr,ZbusOutputStr] = ZbusObj.ReadString(ZbusObj);

% ### Get the models of bus devices
fprintf('Get the descriptor state space model of bus devices.\n')
for i = 1:N_Device
    [GmObj_Cell{i},GmDSS_Cell{i},DevicePara{i},DeviceEqui{i},DeviceDiscreDamping{i},DeviceStateStr{i},DeviceInputStr{i},DeviceOutputStr{i}] = ...
        DeviceModel_Create('Type', DeviceType{i} ,'Flow',PowerFlow{i},'Para',Para{i},'Ts',Ts);
    
    % The following data is not used in the script, but will be used in
    % simulations, please do not delete.
    x_e{i} = DeviceEqui{i}{1};
    u_e{i} = DeviceEqui{i}{2};
    OtherInputs{i} = u_e{i}(3:end,:);
end

% ### Get the model of whole system
fprintf('Get the descriptor state space model of whole system.\n')
GmObj = DeviceModel_Link(GmObj_Cell);
[GsysObj,GsysDSS,Port_v,Port_i,Port_w,Port_T_m,Port_ang_r,Port_P_dc,Port_v_dc] = ...
    GmZbus_Connect(GmObj,ZbusObj);

% ### Chech if the system is proper
fprintf('Check if the whole system is proper:\n')
if isproper(GsysDSS)
    fprintf('Proper.\n');
    fprintf('Calculate the minimum realization of the system model for later use.\n')
    GminSS = minreal(GsysDSS);    
    % This "minreal" function only changes the element sequence of state
    % vectors, but does not change the element sequence of input and output
    % vectors.
    InverseOn = 0;
else
    error('Error: System is improper, which has more zeros than poles.')
end
if is_dss(GminSS)
    error(['Error: Minimum realization is in descriptor state space (dss) form.']);
end

% ### Output the System
fprintf('\n')
fprintf('==================================\n')
fprintf('Print results\n')
fprintf('==================================\n')
fprintf('Model (system object form): GsysObj\n')
fprintf('Model (descriptor state space form): GsysDSS\n')
fprintf('Model (state space form): GminSS\n')
if Enable_PrintOutput
    [SysStateString,SysInputString,SysOutputString] = GsysObj.ReadString(GsysObj);
    fprintf('Print ports of GsysDSS:\n')
    PrintSysString(N_Device,DeviceType,DeviceStateStr,DeviceInputStr,DeviceOutputStr,ZbusStateStr);
	fprintf('Print power flow result:\n')
    fprintf('| bus | P | Q | V | angle | omega |')
    ListPowerFlow
end
    
%%
% ==================================================
% Create Simulink Model
% ==================================================
fprintf('\n')
fprintf('=================================\n')
fprintf('Simulink Model\n')
fprintf('=================================\n')

if Enable_CreateSimulinkModel == 1
    
    fprintf('Creating the simulink model aotumatically, please wait for a few seconds...\n')

    % Set the simulink model name
    Name_Model = 'mymodel_v1';

    % Close existing model with same name
    close_system(Name_Model,0);
    
    % Create the simulink model
    Main_Simulink(Name_Model,ListLine,DeviceType,ListAdvance,PowerFlow);
    fprintf('Get the simulink model successfully.\n')
    fprintf('Warning: for later use of the simulink model, please "save as" a different name.\n')

else
    fprintf('Warning: The auto creation of simulink model is disabled.\n')
end

%%
% ==================================================
% Plot
% ==================================================
    
fprintf('\n')
fprintf('==================================\n')
fprintf('Plot\n')
fprintf('==================================\n')

figure_n = 1000;

% Plot pole/zero map
fprintf('Calculate pole/zero.\n')
pole_sys = pole(GsysDSS)/2/pi;
fprintf('Check if the system is stable:\n')
if isempty(find(real(pole_sys)>1e-9, 1))
    fprintf('Stable.\n');
else
    fprintf('Warning: Unstable.\n')
end
if Enable_PlotPole
    fprintf('Plot pole/zero map.\n')
    figure_n = figure_n+1;
    figure(figure_n);
    scatter(real(pole_sys),imag(pole_sys),'x','LineWidth',1.5); hold on; grid on;
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
    mtit('Global Pole Map');
    
	figure_n = figure_n+1;
    figure(figure_n);
    scatter(real(pole_sys),imag(pole_sys),'x','LineWidth',1.5); hold on; grid on;
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
    mtit('Zoomed Pole Map');
    axis([-250,50,-250,250]);
else
    fprintf('Warning: The default plot of pole map is disabled.\n')
end

omega_p = logspace(-2,4,5e3)*2*pi;
omega_pn = [-flip(omega_p),omega_p];

% Plot admittance
if Enable_PlotAdmittance
    fprintf('Calculate complex-form admittance.\n')
    Tj = [1 1j;     % real to complex transform
          1 -1j];  
    for k = 1:N_Bus
        if DeviceType{k} <= 50
            Gr_ss{k} = GminSS(Port_i([2*k-1,2*k]),Port_v([2*k-1,2*k]));
            Gr_sym{k} = ss2sym(Gr_ss{k});
            Gr_c{k} = Tj*Gr_sym{k}*Tj^(-1);
        end
    end
    fprintf('Plot admittance.\n')
 	figure_n = figure_n+1;
 	figure(figure_n);
    CountLegend = 0;
    VecLegend = {};
    for k = 1:N_Bus
        if DeviceType{k} <= 50
            bodec(Gr_c{k}(1,1),1j*omega_pn,2*pi,'PhaseOn',0); 
            CountLegend = CountLegend + 1;
            VecLegend{CountLegend} = ['Bus',num2str(k)];
        end
    end
    legend(VecLegend);
    mtit('Bode diagram: admittance');
else
    fprintf('Warning: The default plot of admittance spectrum is disabled.\n')
end

% Plot w related
if Enable_PlotSwing
    fprintf('Find the w port relation.\n')
    for k = 1:N_Bus
        if floor(DeviceType{k}/10) == 0
            Gt_ss{k} = GminSS(Port_w(k),Port_T_m(k));
            Gt_sym{k} = -ss2sym(Gt_ss{k});  % The negative sign is because of the load convention.
        elseif floor(DeviceType{k}/10) == 1
         	Gt_ss{k} = GminSS(Port_w(k),Port_ang_r(k));
            Gt_sym{k} = -ss2sym(Gt_ss{k});
        end
    end
    fprintf('Plot w port dynamics.\n')
 	figure_n = figure_n+1;
 	figure(figure_n);
    CountLegend = 0;
    VecLegend = {};
    for k = 1:N_Bus
        if (floor(DeviceType{k}/10) == 0) || (floor(DeviceType{k}/10) == 1)
            bodec(Gt_sym{k},1j*omega_pn,2*pi,'InverseOn',InverseOn,'PhaseOn',0);      
         	CountLegend = CountLegend + 1;
            VecLegend{CountLegend} = ['Bus',num2str(k)]; 
        end
    end
    legend(VecLegend);
    mtit('Bode diagram: omega-related transfer function');
    
else
    fprintf('Warning: The default plot of omega-related spectrum is disabled.\n');
end
    
%%
fprintf('\n')
fprintf('==================================\n')
fprintf('End: toolbox run successfully.\n')
fprintf('==================================\n')
   
%%
% Remained questions:
% - Advanced power flow, including droop bus, etc
% - For system object, please make sure the first port is vdq or idq.
% - The power flow calculation assumes the frequency is Wbase
% - Initialization of network lines, such as the current of line inductor
% and the voltage of line capacitor.
