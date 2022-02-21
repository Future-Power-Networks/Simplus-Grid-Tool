% This function creates a simulink model according to custormer data.

% Author(s): Yitong Li

%% Notes
%
% A very useful function: get_param(gcb,'ObjectParameters')
% The second input argument can be: 'Object Paramters', 'DialogParamters',
% 'PortConnectivity', 'PortHandles', 'ScopeConfiguration'

%%
function MainSimulink(Name_Model,ListBus,ListLine,ApparatusBus,ApparatusType,Advance)

%% Common variables
SimplusGT.Simulink.NewSimulinkModel('ModelName',Name_Model);
% Name_Lib = 'Simplus Grid Tool';
Name_LibFile = 'SimplusGT';
% load_system(Name_LibFile);

%% Shift Base
% The limit when setting position by set_param() is [-32768,32768]. Hence,
% we shift the starting point from [0,0] to a negative coordinate, to give
% more spaces.
Shift_Base = [-30000,-30000];

%% Add powergui
% Parameters
Pos_powergui = [0,0] + Shift_Base;       %
Size_powergui = [70,30];

% Add power gui
SimplusGT.Simulink.SimAddPowerGUI(Name_Model,Size_powergui,Pos_powergui);

%% Add buses
% Parameter
Size_Bus = [30,65];
Pos_Bus{1} = [100,200] + Shift_Base;

% Increase x_distance when the number of mutual branches in creases
[~,MaxCount_ToBus,~] = mode(ListLine(:,2));
Dist_Bus = [200+100*MaxCount_ToBus,300];

% Add bus
[Name_Bus,Pos_Bus] = ...
    SimplusGT.Simulink.SimAddBus(Name_Model,Name_LibFile,Size_Bus,Pos_Bus,ListBus,Dist_Bus);

%% Add apparatuses
% Parameter
Size_Apparatus = [50,90];
Shift_Apparatus = [-150,0];

% Add apparatus
[FullName_Apparatus,Name_Apparatus,Pos_Apparatus] = ...
    SimplusGT.Simulink.SimAddApparatus(Name_Model,Name_LibFile,Size_Apparatus,Shift_Apparatus,Pos_Bus,ApparatusBus,ApparatusType,Advance);

%% Add apparatus ground
% Paramters
Size_D_GND = [20,20];
Shift_D_GND = [20,20];

% Add apparatus ground
Enable_D_GND = 1;
if Enable_D_GND
    SimplusGT.Simulink.SimAddApparatusGround(Name_Model,Size_D_GND,Shift_D_GND,FullName_Apparatus,Name_Apparatus,ApparatusType);
end

%% Add apparatus scope
% Parameters
Size_D_Scope = [30,Size_Apparatus(2)+10];
Shift_D_Scope = [-100,0];
Size_DS_Bus = [5,Size_Apparatus(2)+10];
Shift_DS_Bus = [-30,0];

% Add apparatus scope
SimplusGT.Simulink.SimAddApparatusScope(Name_Model,Size_D_Scope,Shift_D_Scope,Size_DS_Bus,Shift_DS_Bus,Pos_Apparatus,Name_Apparatus,ApparatusType);

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
    SimplusGT.Simulink.SimAddBranch(Name_Model,Name_LibFile,Size_Branch,Shift_Branch,Pos_Bus,ListLine);

% Add transformer
[Name_Trans] = ...
    SimplusGT.Simulink.SimAddTransformer(Name_Model,Size_Trans,Shift_Trans,Pos_Bus,ListLine,Shift_ToBus);

% Add branch ground
SimplusGT.Simulink.SimAddBranchGround(Name_Model,Size_B_GND,Shift_B_GND,FullName_Branch,Name_Branch,ListLine);

% Connect branch to bus
SimplusGT.Simulink.SimConnectBranch2Bus(Name_Model,Name_Bus,Name_Branch,Name_Trans,ListLine);

% Notes:
% Be careful about the settings of ListLine.

%% Connect apparatus to bus
% This procedure is done finally to get a cleaner line auto routing.
SimplusGT.Simulink.SimConnectApparatus2Bus(Name_Model,Name_Bus,Name_Apparatus,ApparatusBus,ApparatusType);

%% Fit the model to view
set_param(gcs,'Zoomfactor','fit to view')

end

