% This function creates a simulink model according to custormer data.

% Author(s): Yitong Li

%% Notes
%
% A very useful function: get_param(gcb,'ObjectParameters')
% The second input argument can be: 'Object Paramters', 'DialogParamters',
% 'PortConnectivity', 'PortHandles', 'ScopeConfiguration'

%%
function MainSimulink(Name_Model,ListLine,DeviceType,ListAdvance,PowerFlow)

%% Common variables
SimplexPS.Simulink.NewSimulinkModel('ModelName',Name_Model);
Name_Lib = 'Simplex Power System';
Name_LibFile = 'SimplexPS';
% load_system(Name_LibFile);

%%
Shift_Base = [-30000,-30000];

%% Add powergui
% Parameters
Pos_powergui = [0,0] + Shift_Base;       %
Size_powergui = [70,30];

% Add power gui
SimplexPS.Simulink.SimAddPowerGUI(Name_Model,Size_powergui,Pos_powergui);

%% Add buses
% Parameter
Size_Bus = [30,65];
Pos_Bus{1} = [100,200] + Shift_Base;

% Increase x_distance when the number of mutual branches in creases
[~,MaxCount_ToBus,~] = mode(ListLine(:,2));
Dist_Bus = [200+100*MaxCount_ToBus,300];

% Add bus
[Name_Bus,Pos_Bus] = ...
    SimplexPS.Simulink.SimAddBus(Name_Model,Name_LibFile,Size_Bus,Pos_Bus,ListLine,Dist_Bus);

%% Add devices
% Parameter
Size_Device = [50,90];
Shift_Device = [-150,0];

% Add device
[FullName_Device,Name_Device,Pos_Device] = ...
    SimplexPS.Simulink.SimAddDevice(Name_Model,Name_LibFile,Size_Device,Shift_Device,Pos_Bus,DeviceType,ListAdvance,PowerFlow);

% Connect device to bus
SimplexPS.Simulink.SimConnectDevice2Bus(Name_Model,Name_Bus,Name_Device,DeviceType);

%% Add device ground
% Paramters
Size_D_GND = [20,20];
Shift_D_GND = [20,20];

% Add device ground
Enable_D_GND = 1;
if Enable_D_GND
    SimplexPS.Simulink.SimAddDeviceGround(Name_Model,Size_D_GND,Shift_D_GND,FullName_Device,Name_Device,DeviceType);
end

%% Add device scope
% Parameters
Size_D_Scope = [30,Size_Device(2)+10];
Shift_D_Scope = [-100,0];
Size_DS_Bus = [5,Size_Device(2)+10];
Shift_DS_Bus = [-30,0];

% Add device scope
SimplexPS.Simulink.SimAddDeviceScope(Name_Model,Size_D_Scope,Shift_D_Scope,Size_DS_Bus,Shift_DS_Bus,Pos_Device,Name_Device,DeviceType);

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
    SimplexPS.Simulink.SimAddBranch(Name_Model,Size_Branch,Shift_Branch,Pos_Bus,ListLine);

% Add transformer
[Name_Trans] = ...
    SimplexPS.Simulink.SimAddTransformer(Name_Model,Size_Trans,Shift_Trans,Pos_Bus,ListLine,Shift_ToBus);

% Add branch ground
SimplexPS.Simulink.SimAddBranchGround(Name_Model,Size_B_GND,Shift_B_GND,FullName_Branch,Name_Branch,ListLine);

% Connect branch to bus
SimplexPS.Simulink.SimConnectBranch2Bus(Name_Model,Name_Bus,Name_Branch,Name_Trans,ListLine);

%% Fit the model to view
% set_param(gcs, 'ZoomFactor','FitSystem')
set_param(gcs,'Zoomfactor','fit to view')

end

