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
% - The connection of inductor and SG, or the port type of SG. A simple
% method: provide two types of port: v -> i and i-> v
% - Put the base value into excel
% - Maybe it is better to put default parameters into excel
% - Deal with the transformer
% - Bus, device date, should be re-arranged somewhere, so that they can be
% in order

% The device from matlab is also in current-source type or voltage-source
% type.

%%
% ==================================================
% Common paramters
% ==================================================
fprintf('Get common parameters.\n')
% Sampling for simulation
Fs = 50e3;          % (Hz) samping frequency
Ts = 1/Fs;          % (s) Samping period

% Angular frequency
W0 = 2*pi*50;     	% (rad/s)

Vbase = 1;
Sbase = 1;
Ibase = Sbase/Vbase;
Wbase = W0;

%%
% ==================================================
% Load customized data
% ==================================================
fprintf('Load customized data.\n')

% ### Load the data
ListBus = xlsread('netlist.xlsx',1);     % Other possible function: readmatrix, csvread ...
ListLine = xlsread('netlist.xlsx',2);
ListDevice = xlsread('netlist.xlsx',3);

% ### Re-arrange the netlist and check error
[ListBus,ListLine,ListDevice,N_Bus,N_Branch,N_Device] = RearrangeNetlist(ListBus,ListLine,ListDevice);

% ### Re-arrange the device data
[DeviceType,Para] = RearrangeDeviceData(ListDevice,W0);

%%
% ==================================================
% Descriptor state space model
% ==================================================

% ### Power flow analysis
fprintf('Do the power flow analysis.\n')
[PowerFlow,~,~,~,~,~,~]=PowerFlow_GS(ListBus,ListLine,W0);

% ### Get the model of lines
fprintf('Get the descriptor-state-space model of network lines.\n')
[YbusObj,~,~] = YbusCalcDSS(ListLine,W0);
[~,YbusDSS] = YbusObj.ReadDSS(YbusObj);
[~,lsw] = size(YbusDSS.B);
ZbusObj = obj_SwitchInOut(YbusObj,lsw);

% ### Get the models of bus devices
fprintf('Get the descriptor-state-space model of bus devices.\n')
for i = 1:N_Device
    [GmObj_Cell{i},GmDSS_Cell{i},DevicePara{i},DeviceEqui{i}] = ...
        DeviceModel_Create('Type', DeviceType{i} ,'Flow',PowerFlow{i},'Para',Para{i});
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
    G_SS = minreal(GsysDSS);    
    % This function only changes the element sequence of state vectors, but
    % does not change the element sequence of input and output vectors.
    InverseOn = 0;
elseif isproper(1/GsysDSS)
    fprintf('Warrning: improper.\n');
    fprintf('Calculate the minimum realization of the inverse of the system model for later use.\n')
    G_SS = minreal(1/GsysDSS);
    InverseOn = 1;
else
    error(['Error: both the system and the its inverse are improper.'])
end
if ~isempty(G_SS.E)
    error(['Error: minimum realization is in descriptor state space form.']);
end

% ### Output the System
fprintf('### Output the system\n')
fprintf('System obj name: GsysObj\n')
[StateString,InputStr,OutputStr] = GsysObj.ReadString(GsysObj)
fprintf('Minimum realization system name: G_SS\n')

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
main_simulink(Name_Model,W0,ListLine,N_Bus,N_Branch,N_Device,DeviceType);
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
            Gr_ss{k} = G_SS(iport([2*k-1,2*k]),vport([2*k-1,2*k]));
        else          
            Gr_ss{k} = G_SS(vport([2*k-1,2*k]),iport([2*k-1,2*k]));
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
   
