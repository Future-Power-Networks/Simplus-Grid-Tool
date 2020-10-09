% Author(s): Yitong Li, Yunjie Gu

clear all; 
clc;
close all;

% cmap = get(groot,'defaultAxesColorOrder');
fprintf('### Start.\n')

%%
% check1 = 0      % For test, will be deleted later
% Remained questions:
% - Speed up the power flow
% - Advanced power flow, including droop bus, etc
% - Put the base value into excel

%%
% ==================================================
% Load customized data
% ==================================================
fprintf('Load customized data.\n')

% ### Load the data
% Other possible function: readmatrix, csvread ...
Name_Netlist = 'netlist.xlsx';
% Name_Netlist = 'netlist_TestToolbox.xlsx';
% Name_Netlist = 'netlist_TestInfiniteBus.xlsx';
ListBus = xlsread(Name_Netlist,1);     
ListDevice = xlsread(Name_Netlist,2);
ListSimulation = xlsread(Name_Netlist,3);
ListLine = xlsread(Name_Netlist,4);
ListLineIEEE = xlsread(Name_Netlist,5);

% ### Re-arrange the simulation data
Fs = ListSimulation(1);
Ts = 1/Fs;  % (s) sampling period
FundamentalFreq = ListSimulation(length(ListSimulation));
W0 = FundamentalFreq*2*pi;

% ### Re-arrange the netlist and check error
[ListLine] = RearrangeNetlist_IEEE2Toolbox(ListLine,ListLineIEEE);
[ListBus,ListLine,ListDevice,N_Bus,N_Branch,N_Device] = RearrangeNetlist(ListBus,ListLine,ListDevice);

% ### Re-arrange the device data
[DeviceType,Para] = RearrangeDeviceData(ListDevice,W0);

%%
% ==================================================
% Common paramters
% ==================================================
fprintf('Get common parameters.\n')
% Base values
Vbase = 1;
Sbase = 1;
Ibase = Sbase/Vbase;
Wbase = W0;

%%
% ==================================================
% Descriptor state space model
% ==================================================

% ### Power flow analysis
fprintf('Do the power flow analysis.\n')
[PowerFlow,~,~,~,~,~,~]=PowerFlow_GS(ListBus,ListLine,W0);

% ### Get the model of lines
fprintf('Get the descriptor-state-space model of network lines.\n')
% Move load flow to bus admittance matrix
[ListBus,ListLine,PowerFlow] = Load2SelfBranch(ListBus,ListLine,DeviceType,PowerFlow);
[YbusObj,~,~] = YbusCalcDSS(ListLine,W0);
[~,YbusDSS] = YbusObj.ReadDSS(YbusObj);
[~,lsw] = size(YbusDSS.B);
ZbusObj = obj_SwitchInOut(YbusObj,lsw);

% ### Get the models of bus devices
fprintf('Get the descriptor-state-space model of bus devices.\n')
for i = 1:N_Device
    [GmObj_Cell{i},GmDSS_Cell{i},DevicePara{i},DeviceEqui{i},DeviceDiscreDamping{i}] = ...
        DeviceModel_Create('Type', DeviceType{i} ,'Flow',PowerFlow{i},'Para',Para{i},'Ts',Ts);
    
    % The following data is not used in the script, but will be used in
    % simulations, please do not delete.
    x_e{i} = DeviceEqui{i}{1};
    u_e{i} = DeviceEqui{i}{2};
    OtherInputs{i} = u_e{i}(3:end,:);
end

% ### Get the model of whole system
fprintf('Get the descriptor-state-space model of whole system.\n')
GmObj = DeviceModel_Link(GmObj_Cell);
[GsysObj,GsysDSS,vport,iport] = GmZbus_Connect(GmObj,ZbusObj);

% ### Chech if the system is proper
fprintf('Check if the whole system is proper:\n')
if isproper(GsysDSS)
    fprintf('Proper.\n');
    fprintf('Calculate the minimum realization of the system model for later use.\n')
    GminSS = minreal(GsysDSS);    
    % This function only changes the element sequence of state vectors, but
    % does not change the element sequence of input and output vectors.
    InverseOn = 0;
elseif isproper(1/GsysDSS)
    fprintf('Warrning: improper.\n');
    fprintf('Calculate the minimum realization of the inverse of the system model for later use.\n')
    GminSS = minreal(1/GsysDSS);
    InverseOn = 1;
else
    error(['Error: both the system and the its inverse are improper.'])
end
if ~isempty(GminSS.E)
    error(['Error: minimum realization is in descriptor state space form.']);
end

% ### Output the System
fprintf('### Output the system\n')
fprintf('System object name: GsysObj\n')
fprintf('System name: GsysDSS\n')
fprintf('Minimum realization system name: GminSS\n')
[StateString,InputString,OutputString] = GsysObj.ReadString(GsysObj)

%%
% ==================================================
% Create Simulink Model
% ==================================================
fprintf('### Simulink Model\n')
fprintf('Create the simulink model aotumatically.\n')

% Set the simulink model name
Name_Model = 'mymodel_v1';

% Close existing model with same name
close_system(Name_Model,0);

% Create the simulink model
Main_Simulink(Name_Model,ListLine,DeviceType,ListSimulation,PowerFlow);
fprintf('Get the simulink model successfully.\n')
fprintf('Warning: for later use of the simulink model, please "save as" a different name.\n')

%%
% ==================================================
% Plot
% ==================================================
    
fprintf('### Plot\n')

figure_n = 0;

% Plot pole/zero map
if 1
    figure_n = figure_n+1;
    figure(figure_n);
    fprintf('Calculate pole/zero.\n')
    psys = pole(GsysDSS)/2/pi;
    fprintf('Plot pole/zero map.\n')
    scatter(real(psys),imag(psys),'x','LineWidth',1.5); hold on; grid on;
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
    axis([-25,5,-80,80]);
end

omega_p = logspace(-2,3,1e3)*2*pi;
omega_pn = [-flip(omega_p),omega_p];

% Plot admittance: method 1
if 1 
    fprintf('Calculate admittance.\n')
    Tj = [1 1j;     % real to complex transform
          1 -1j];  
    for k = 1:N_Bus
        if InverseOn == 0
            Gr_ss{k} = GminSS(iport([2*k-1,2*k]),vport([2*k-1,2*k]));
        else          
            Gr_ss{k} = GminSS(vport([2*k-1,2*k]),iport([2*k-1,2*k]));
        end
        Gr_sym{k} = ss2sym(Gr_ss{k});
        Gr_c{k} = Tj*Gr_sym{k}*Tj^(-1);
    end
    fprintf('Plot admittance.\n')
 	figure_n = figure_n+1;
 	figure(figure_n);
    for k = 1:N_Bus
        bodec(Gr_c{k}(1,1),1j*omega_pn,2*pi,'InverseOn',InverseOn,'PhaseOn',0);                           
    end
end
    
%%
fprintf('### End: toolbox run successfully.\n')
   
