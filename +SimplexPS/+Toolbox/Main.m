% Author(s): Yitong Li, Yunjie Gu

%% 
% Notes:
%
% Please read "README.md" first before using the toolbox.
%
% Please use "Main_Customer.m" rather than this file for running the
% toolbox.

%%
fprintf('==================================\n')
fprintf('Start to run SimplexPowerSystem\n')
fprintf('==================================\n')

%%
% ==================================================
% Load customized data
% ==================================================
fprintf('Loading data from "UserData.xlsx", please wait a second...\n')

% ### Re-arrange basic settings
ListBasic    = xlsread(UserData,'Basic');
Fs = ListBasic(1);
Ts = 1/Fs;              % (s), sampling period
Fbase = ListBasic(2);   % (Hz), base frequency
Sbase = ListBasic(3);   % (VA), base power
Vbase = ListBasic(4);   % (V), base voltage
Ibase = Sbase/Vbase;    % (A), base current
Zbase = Sbase/Ibase;    % (Ohm), base impedance
Ybase = 1/Zbase;        % (S), base admittance
Wbase = Fbase*2*pi;     % (rad/s), base angular frequency

% ### Re-arrange advanced settings
ListAdvance  = xlsread(UserData,'Advance');
Flag_PowerFlowAlgorithm      	= ListAdvance(5);
Enable_CreateSimulinkModel      = ListAdvance(6);
Enable_PlotPole                 = ListAdvance(7);
Enable_PlotAdmittance           = ListAdvance(8);
Enable_PrintOutput              = ListAdvance(9);
Enable_PlotSwing                = ListAdvance(10);
if length(ListAdvance) ==10
    Enable_Participation = 0;
else
    Enable_Participation = ListAdvance(11);
end

% ### Re-arrange the bus netlist
[ListBus,N_Bus] = SimplexPS.Toolbox.RearrangeListBus(UserData);

% ### Re-arrange the line netlist
[ListLine,N_Branch.N_Bus_] = SimplexPS.Toolbox.RearrangeListLine(UserData,ListBus);

% ### Re-arrange the device netlist
[DeviceBus,DeviceType,Para,N_Device] = SimplexPS.Toolbox.RearrangeListDevice(UserData,Wbase,ListBus);
% The names of "DeviceType" and "Para" can not be changed, because they
% will also be used in simulink model.

% Notes:
% No error checking if number of devices is different from number of buses.

%%
% ==================================================
% Power flow analysis
% ==================================================

% ### Power flow analysis
fprintf('Doing the power flow analysis...\n')
switch Flag_PowerFlowAlgorithm
    case 1  % Gauss-Seidel 
        [PowerFlow,~,~,~,~,~,~,~] = SimplexPS.PowerFlow.PowerFlowGS(ListBus,ListLine,Wbase);
    case 2  % Newton-Raphson
        [PowerFlow] = SimplexPS.PowerFlow.PowerFlowNR(ListBus,ListLine,Wbase);
    otherwise
        error(['Error: Wrong setting for power flow algorithm.']);
end
% Form of PowerFlow{i}: P, Q, V, xi, w

% For printing later
ListPowerFlow = SimplexPS.PowerFlow.Rearrange(PowerFlow);

% Move load flow (PLi and QLi) to bus admittance matrix
[ListBus,ListLineNew,PowerFlowNew] = SimplexPS.PowerFlow.Load2SelfBranch(ListBus,ListLine,PowerFlow);

% For printting later
ListPowerFlow_ = SimplexPS.PowerFlow.Rearrange(PowerFlowNew);

%%
% ==================================================
% Descriptor state space model
% ==================================================

% ### Get the model of lines
fprintf('Getting the descriptor state space model of network lines...\n')

[YbusObj,YbusDSS,~] = SimplexPS.Toolbox.YbusCalcDss(ListBus,ListLineNew,Wbase);
[~,lsw] = size(YbusDSS.B);
ZbusObj = SimplexPS.ObjSwitchInOut(YbusObj,lsw);
[ZbusStateStr,ZbusInputStr,ZbusOutputStr] = ZbusObj.GetString(ZbusObj);

% ### Get the models of bus devices
fprintf('Getting the descriptor state space model of bus devices...\n')
for i = 1:N_Device
    if length(DeviceBus{i}) == 1
     	DevicePowerFlow{i} = PowerFlowNew{DeviceBus{i}};
    elseif length(DeviceBus{i}) == 2
        DevicePowerFlow{i} = [PowerFlowNew{DeviceBus{i}(1)},PowerFlowNew{DeviceBus{i}(2)}];
    else
        error(['Error']);
    end
    [GmObj_Cell{i},GmDSS_Cell{i},DevicePara{i},DeviceEqui{i},DeviceDiscreDamping{i},OtherInputs{i}] = ...
        SimplexPS.Toolbox.DeviceModelCreate(DeviceBus{i},DeviceType{i},DevicePowerFlow{i},Para{i},Ts,ListBus);
    
    % The following data is not used in the script, but will be used in
    % simulations. Do not delete!
    x_e{i} = DeviceEqui{i}{1};
    u_e{i} = DeviceEqui{i}{2};
end

% ### Get the appended model of all devices
fprintf('Getting the appended descriptor state space model of all devices...\n')
GmObj = SimplexPS.Toolbox.DeviceModelLink(GmObj_Cell);

% ### Get the model of whole system
fprintf('Getting the descriptor state space model of whole system...\n')
[GsysObj,GsysDSS,Port_v,Port_i,Port_w,Port_T_m,Port_ang_r,Port_P_dc,Port_v_dc] = ...
    SimplexPS.Toolbox.ConnectGmZbus(GmObj,ZbusObj,N_Bus);

% ### Whole-system admittance model
YsysDSS = GsysDSS(Port_i,Port_v);   

% ### Chech if the system is proper
fprintf('Checking if the whole system is proper:\n')
if isproper(GsysDSS)
    fprintf('Proper!\n');
    fprintf('Calculating the minimum realization of the system model for later use...\n')
    % GminSS = minreal(GsysDSS);
    GminSS = SimplexPS.dss2ss(GsysDSS);
    % This "minreal" function only changes the element sequence of state
    % vectors, but does not change the element sequence of input and output
    % vectors.
    InverseOn = 0;
else
    error('Error: System is improper, which has more zeros than poles.')
end
if SimplexPS.is_dss(GminSS)
    error(['Error: Minimum realization is in descriptor state space (dss) form.']);
end

%%
% ==================================================
% Print results
% ==================================================

% ### Output the System
fprintf('\n')
fprintf('==================================\n')
fprintf('Print results\n')
fprintf('==================================\n')
fprintf('Model (system object form): GsysObj\n')
fprintf('Model (descriptor state space form): GsysDSS\n')
if Enable_PrintOutput
    [SysStateString,SysInputString,SysOutputString] = GsysObj.GetString(GsysObj);
    fprintf('Print ports of GsysDSS:\n')
    SimplexPS.Toolbox.PrintSysString(DeviceBus,DeviceType,GmObj_Cell,ZbusStateStr);
	fprintf('Print power flow result:\n')
    fprintf('| bus | P | Q | V | angle | omega |')
    ListPowerFlow
end

fprintf('Other models: \n')
fprintf('Model (state space form): GminSS\n')
fprintf('Model (descriptor state space form): YsysDSS\n')

%%
% ==================================================
% Check Stability
% ==================================================

fprintf('\n')
fprintf('==================================\n')
fprintf('Check Stability\n')
fprintf('==================================\n')

fprintf('Calculatting pole/zero...\n')
pole_sys = pole(GsysDSS)/2/pi;
fprintf('Checking if the system is stable:\n')
if isempty(find(real(pole_sys)>1e-9, 1))
    fprintf('Stable!\n');
else
    fprintf('Warning: Unstable!\n')
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
    
    fprintf('Creating the simulink model automatically, please wait a second...\n')

    % Set the simulink model name
    Name_Model = 'mymodel_v1';

    % Close existing model with same name
    close_system(Name_Model,0);
    
    % Create the simulink model
    SimplexPS.Simulink.MainSimulink(Name_Model,ListBus,ListLineNew,DeviceBus,DeviceType,ListAdvance,PowerFlowNew);
    fprintf('Get the simulink model successfully! \n')
    fprintf('Please click the "run" button in the model to run it.\n')
    %fprintf('Warning: for later use of the simulink model, please "save as" a different name.\n')

else
    fprintf('Warning: The auto creation of simulink model is disabled.\n')
end

%%
% ==================================================
% Plot
% ==================================================

fprintf('\n')
fprintf('==================================\n')
fprintf('Plot Fundamentals\n')
fprintf('==================================\n')

% Initialize figure index
figure_n = 1000;

% Plot pole/zero map
if Enable_PlotPole
    fprintf('Plotting pole map...\n')
    figure_n = figure_n+1;
    figure(figure_n);
    subplot(1,2,1)
    scatter(real(pole_sys),imag(pole_sys),'x','LineWidth',1.5); hold on; grid on;
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
    title('Global pole map');
    
	subplot(1,2,2)
    scatter(real(pole_sys),imag(pole_sys),'x','LineWidth',1.5); hold on; grid on;
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
    title('Zoomed pole map');
    axis([-80,20,-150,150]);
    
    %SimplexPS.mtit('Pole Map');
else
    fprintf('Warning: The default plot of pole map is disabled.\n')
end

stop

% ====================================================================
% The following plot functions have to be re-written!
% ====================================================================

omega_p = logspace(-2,4,5e3)*2*pi;
omega_pn = [-flip(omega_p),omega_p];

% Plot admittance
if Enable_PlotAdmittance
    fprintf('Plotting admittance spectrum...\n')
    Tj = [1 1j;     % real to complex transform
          1 -1j];  
    for k = 1:N_Bus
        if DeviceType{k} <= 50
            Gr_ss{k} = GminSS(Port_i([2*k-1,2*k]),Port_v([2*k-1,2*k]));
            Gr_sym{k} = SimplexPS.ss2sym(Gr_ss{k});
            Gr_c{k} = Tj*Gr_sym{k}*Tj^(-1);
        end
    end
 	figure_n = figure_n+1;
 	figure(figure_n);
    CountLegend = 0;
    VecLegend = {};
    for k = 1:N_Bus
        if DeviceType{k} <= 50
            SimplexPS.bode_c(Gr_c{k}(1,1),1j*omega_pn,2*pi,'PhaseOn',0); 
            CountLegend = CountLegend + 1;
            VecLegend{CountLegend} = ['Bus',num2str(k)];
        end
    end
    legend(VecLegend);
    SimplexPS.mtit('Bode diagram: admittance');
else
    fprintf('Warning: The default plot of admittance spectrum is disabled.\n')
end

% Plot w related
if Enable_PlotSwing
    fprintf('Plotting frequency-port dynamics...\n')
    for k = 1:N_Bus
        if floor(DeviceType{k}/10) == 0
            Gt_ss{k} = GminSS(Port_w(k),Port_T_m(k));
            Gt_sym{k} = -SimplexPS.ss2sym(Gt_ss{k});  % The negative sign is because of the load convention.
        elseif floor(DeviceType{k}/10) == 1
         	Gt_ss{k} = GminSS(Port_w(k),Port_ang_r(k));
            Gt_sym{k} = -SimplexPS.ss2sym(Gt_ss{k});
        end
    end
 	figure_n = figure_n+1;
 	figure(figure_n);
    CountLegend = 0;
    VecLegend = {};
    for k = 1:N_Bus
        if (floor(DeviceType{k}/10) == 0) || (floor(DeviceType{k}/10) == 1)
            SimplexPS.bode_c(Gt_sym{k},1j*omega_pn,2*pi,'PhaseOn',0);      
         	CountLegend = CountLegend + 1;
            VecLegend{CountLegend} = ['Bus',num2str(k)]; 
        end
    end
    legend(VecLegend);
    SimplexPS.mtit('Bode diagram: frequency-port transfer function');
    
else
    fprintf('Warning: The default plot of omega-related spectrum is disabled.\n');
end

%%
fprintf('\n')
fprintf('==================================\n')
fprintf('Modal Analysis\n')
fprintf('==================================\n')
if Enable_Participation == 1
    SimplexPS.Modal.ModalPreRun;
    SimplexPS.Modal.ModalAnalysis;
    fprintf('Generating GreyboxConfg.xlsx for user to config Greybox analysis.\n');    
else
    fprintf('Warning: The modal (participation) analysis is disabled.');
end


%%
fprintf('\n')
fprintf('==================================\n')
fprintf('End: run successfully.\n')
fprintf('==================================\n')
%%
% Remained questions:
% - Advanced power flow, including droop bus, etc
% - For system object, please make sure the first port is vdq or idq.
% - The power flow calculation assumes the frequency is Wbase
% - Initialization of network lines, such as the current of line inductor
% and the voltage of line capacitor.
