% ### Re-arrange basic settings
ListBasic = xlsread(UserData,'Basic');

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
ListAdvance = xlsread(UserData,'Advance');
Flag_PowerFlowAlgorithm   	= ListAdvance(5);
Enable_CreateSimulinkModel	= ListAdvance(6);
Enable_PlotPole           	= ListAdvance(7);
Enable_PlotAdmittance     	= ListAdvance(8);
Enable_PrintOutput       	= ListAdvance(9);
Enable_Participation        = ListAdvance(10);

% ### Re-arrange the bus netlist
[ListBus,N_Bus] = SimplexPS.Toolbox.RearrangeListBus(UserData);

% ### Re-arrange the line netlist
[ListLine,N_Branch.N_Bus_] = SimplexPS.Toolbox.RearrangeListLine(UserData,ListBus);

% ### Re-arrange the device netlist
[DeviceBus,DeviceType,DevicePara,N_Device] = SimplexPS.Toolbox.RearrangeListDevice(UserData,Wbase,ListBus);

% Notes:
% The codes in this part are borrowed from the SimplexPS toolbox.