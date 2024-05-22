UserData = 'IEEE_39Bus.xlsx';
% Change current path of matlab
% ==================================================
% pathstr = mfilename('fullpath');        % Get the path of main.m
% [pathstr,~,~]  = fileparts(pathstr);    % Get the path of Toolbox namespace
% [pathstr,~,~]  = fileparts(pathstr);    % Get the of root namespace
% [pathstr,~,~]  = fileparts(pathstr);    % Get the path of toolbox
% cd(pathstr);                            % Change the current address

%%
% ==================================================
% Load customized data
% ==================================================
% fprintf('Loading data, please wait a second...\n')

% ### Re-arrange basic settings
ListBasic = xlsread(UserData,'Basic');
Fs = ListBasic(1);
Ts = 1/Fs;              % (s), sampling period
Fbase = ListBasic(2);   % (Hz), base frequency
Sbase = ListBasic(3);   % (VA), base power
Vbase = ListBasic(4);   % (V), base voltage
Ibase = Sbase/Vbase;    % (A), base current
Zbase = Vbase/Ibase;    % (Ohm), base impedance
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
[ListBus,N_Bus] = SimplusGT.Toolbox.RearrangeListBus(UserData);

ListBus_=ListBus; % temparary added by Yue to save the original data

% ### Re-arrange the line netlist
[ListLine,N_Branch.N_Bus_] = SimplusGT.Toolbox.RearrangeListLine(UserData,ListBus);
DcAreaFlag = find(ListBus(:,12)==2);

% ### Re-arrange the apparatus netlist
[ApparatusBus,ApparatusType,Para,N_Apparatus] = SimplusGT.Toolbox.RearrangeListApparatus(UserData,Wbase,ListBus);
% The names of "ApparatusType" and "Para" can not be changed, because they
% will also be used in simulink model.

% Notes:
% No error checking if number of apparatuses is different from number of buses.

%%
% ==================================================
% Power flow analysis
% ==================================================

% ### Power flow analysis
fprintf('Doing the power flow analysis...\n')
if ~isempty(DcAreaFlag)
    Flag_PowerFlowAlgorithm = 1;
    fprintf(['Warning: Because the system has dc area(s), the Gauss-Seidel power flow method is always used.\n']);
end
switch Flag_PowerFlowAlgorithm
    case 1  % Gauss-Seidel 
        [PowerFlow,~,~,~,~,~,~,~] = SimplusGT.PowerFlow.PowerFlowGS(ListBus,ListLine,Wbase);
    case 2  % Newton-Raphson
       	[PowerFlow] = SimplusGT.PowerFlow.PowerFlowNR(ListBus,ListLine,Wbase);
    otherwise
        error(['Error: Wrong setting for power flow algorithm.']);
end
% Form of PowerFlow{i}: P, Q, V, xi, w
% P and Q are in load convention, i.e., the P and Q flowing from the bus to
% the apparatus.

% For printing later
ListPowerFlow = SimplusGT.PowerFlow.Rearrange(PowerFlow);

YN = SimplusGT.PowerFlow.YbusCalc(ListLine);
YA = zeros(39);
YA(1,1) = 1/(1i*0.296);
YN_ = YN+YA;
Zth_matrix = inv(YN_);
Zth_30_38 = abs(diag(Zth_matrix(30:38,30:38)));
Zth_30_38
%YA()

PQ_load = ListBus(:,7:8);         % PLi, QLi..
Y_load = diag(PQ_load(:,1)-1i*PQ_load(:,2));
YN_all = YN+YA+Y_load;
Zth_matrix_all = inv(YN_all);
Zth_30_38_all = abs(diag(Zth_matrix_all(30:38,30:38)));
Zth_30_38_all


