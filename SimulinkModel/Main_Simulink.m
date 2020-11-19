% This function creates a simulink model according to custormer data.

% Author(s): Yitong Li

%% Notes
%
% A very useful function: get_param(gcb,'ObjectParameters')
% The second input argument can be: 'Object Paramters', 'DialogParamters',
% 'PortConnectivity', 'PortHandles', 'ScopeConfiguration'

%%
function Main_Simulink(Name_Model,ListLine,DeviceType,ListAdvance,PowerFlow)

%% Common variables
NewSimulinkModel('ModelName',Name_Model);
Name_Lib = 'Simplex Power Systems';
Name_LibFile = 'SimplexPS_2016a';
% load_system(Name_LibFile);

%% Add powergui
% Parameters
Pos_powergui = [0,0];       %
Size_powergui = [70,30];

% Add power gui
Sim_AddPowerGUI(Name_Model,Size_powergui,Pos_powergui);

%% Add buses
% Parameter
Size_Bus = [30,65];
Pos_Bus{1} = [100,200];

% Increase x_distance when the number of mutual branches in creases
[~,MaxCount_ToBus,~] = mode(ListLine(:,2));
Dist_Bus = [200+100*MaxCount_ToBus,300];

% Add bus
[Name_Bus,Pos_Bus] = ...
    Sim_AddBus(Name_Model,Name_LibFile,Size_Bus,Pos_Bus,ListLine,Dist_Bus);

%% Add devices
% Parameter
Size_Device = [50,90];
Shift_Device = [-150,0];

% Add device
[FullName_Device,Name_Device,Pos_Device] = ...
    Sim_AddDevice(Name_Model,Name_LibFile,Size_Device,Shift_Device,Pos_Bus,DeviceType,ListAdvance,PowerFlow);

% Connect device to bus
Sim_ConnectDevice2Bus(Name_Model,Name_Bus,Name_Device,DeviceType);

%% Add device ground
% Paramters
Size_D_GND = [20,20];
Shift_D_GND = [20,20];

% Add device ground
Enable_D_GND = 1;
if Enable_D_GND
    Sim_AddDeviceGround(Name_Model,Size_D_GND,Shift_D_GND,FullName_Device,Name_Device,DeviceType);
end

%% Add device scope
% Parameters
Size_D_Scope = [30,Size_Device(2)+10];
Shift_D_Scope = [-100,0];
Size_DS_Bus = [5,Size_Device(2)+10];
Shift_DS_Bus = [-30,0];

% Add device scope
Sim_AddDeviceScope(Name_Model,Size_D_Scope,Shift_D_Scope,Size_DS_Bus,Shift_DS_Bus,Pos_Device,Name_Device,DeviceType);

%% Add branches
% Parameter
Size_Branch = [Size_Bus(2),Size_Bus(2)];    % Branch
Shift_Branch  = [+100,+100];
Size_Trans = Size_Branch;                   % Transformer
Shift_Trans = Shift_Branch;
Size_B_GND = Size_D_GND;                    % Ground
Shift_B_GND = [-Size_D_GND(1)/2,30];

% Add branch
[FullName_Branch,Name_Branch,Shift_ToBus] = ...
    Sim_AddBranch(Name_Model,Size_Branch,Shift_Branch,Pos_Bus,ListLine);

% Add transformer
[Name_Trans] = ...
    Sim_AddTransformer(Name_Model,Size_Trans,Shift_Trans,Pos_Bus,ListLine,Shift_ToBus);

% Add branch ground
Sim_AddBranchGround(Name_Model,Size_B_GND,Shift_B_GND,FullName_Branch,Name_Branch,ListLine);

% Connect branch to bus
Sim_ConnectBranch2Bus(Name_Model,Name_Bus,Name_Branch,Name_Trans,ListLine);

end

