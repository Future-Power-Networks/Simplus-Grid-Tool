%
% This script sets up a model in preparation for running
% This uses a Matlab struct to define the model parameters (inputData)
% which needs to be loaded from a json file
%
% Author(s): Rob Oldaker
%

%%
fprintf('==================================\n')
fprintf('Start to run Simplus Grid Tool\n')
fprintf('==================================\n')

%% 
% ==================================================
% Change current path of matlab
% ==================================================
pathstr = mfilename('fullpath');        % Get the path of main.m
[pathstr,~,~]  = fileparts(pathstr);    % Get the path of Toolbox namespace
[pathstr,~,~]  = fileparts(pathstr);    % Get the of root namespace
[pathstr,~,~]  = fileparts(pathstr);    % Get the path of toolbox
cd(pathstr);                            % Change the current address

%%
% ==================================================
% Load customized data
% ==================================================
fprintf('Loading data, please wait a second...\n')

% ### Re-arrange basic settings
Basic = inputData.basic;
Fs = Basic.Fs;
Ts = 1/Fs;               % (s), sampling period
Fbase = Basic.Fbase; % (Hz), base frequency
Sbase = Basic.Sbase; % (VA), base power
Vbase = Basic.Vbase; % (V), base voltage
Ibase = Sbase/Vbase;     % (A), base current
Zbase = Vbase/Ibase;     % (Ohm), base impedance
Ybase = 1/Zbase;         % (S), base admittance
Wbase = Fbase*2*pi;      % (rad/s), base angular frequency

% ### Re-arrange advanced settings
advData = inputData.adv;
Flag_PowerFlowAlgorithm   	= advData.powerFlowAlgorithm;
Enable_CreateSimulinkModel	= advData.enableCreateSimulinkModel;
Enable_PlotPole           	= advData.enablePlotPole;
Enable_PlotAdmittance     	= advData.enablePlotAdmittance;
Enable_PrintOutput       	= advData.enablePrintOutput;
Enable_Participation        = advData.enableParticipation;
ListAdvance = [];
ListAdvance(1)=advData.discretizationMethod;
ListAdvance(2)=advData.linearizationTimes;
ListAdvance(3)=advData.discretizationDampingFlag;
ListAdvance(4)=advData.directFeedThrough;
ListAdvance(5)=advData.powerFlowAlgorithm;
ListAdvance(6)=advData.enableCreateSimulinkModel;
ListAdvance(7)=advData.enablePlotPole;
ListAdvance(8)=advData.enablePlotAdmittance;
ListAdvance(9)=advData.enablePrintOutput;
ListAdvance(10)=advData.enableParticipation;
ListAdvance=ListAdvance';

% ### Re-arrange the bus netlist
[ListBus,N_Bus] = SimplusGT.Toolbox.RearrangeListBusStruct(inputData);

% ### Re-arrange the line netlist
[ListLine,N_Branch.N_Bus_] = SimplusGT.Toolbox.RearrangeListLineStruct(inputData,ListBus);
DcAreaFlag = find(ListBus(:,12)==2);

% ### Re-arrange the apparatus netlist
[ApparatusBus,ApparatusType,Para,N_Apparatus] = SimplusGT.Toolbox.RearrangeListApparatusStruct(inputData,Wbase,ListBus);
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
        [PowerFlow] = SimplusGT.PowerFlow.PowerFlowGS(ListBus,ListLine,Wbase);
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

% Move load flow (PLi and QLi) to bus admittance matrix
[ListBusNew,ListLineNew,PowerFlowNew] = SimplusGT.PowerFlow.Load2SelfBranch(ListBus,ListLine,PowerFlow);

% For printting later
ListPowerFlowNew = SimplusGT.PowerFlow.Rearrange(PowerFlowNew);

%%
% ==================================================
% Descriptor state space model
% ==================================================

% ### Get the model of lines
fprintf('Getting the descriptor state space model of network lines...\n')

[YbusObj,YbusDSS,~] = SimplusGT.Toolbox.YbusCalcDss(ListBusNew,ListLineNew,Wbase);
% if Enable_RemoveDuplicateStates == 1
%     YbusObj = SimplusGT.Toolbox.RemoveYbusDuplicateStates(YbusObj);
% end
ZbusObj = SimplusGT.ObjSwitchInOut(YbusObj,length(YbusDSS));

% ### Get the models of bus apparatuses
fprintf('Getting the descriptor state space model of bus apparatuses...\n')
for i = 1:N_Apparatus
    if length(ApparatusBus{i}) == 1
     	ApparatusPowerFlow{i} = PowerFlowNew{ApparatusBus{i}};
    elseif length(ApparatusBus{i}) == 2
        ApparatusPowerFlow{i} = [PowerFlowNew{ApparatusBus{i}(1)},PowerFlowNew{ApparatusBus{i}(2)}];
    else
        error(['Error']);
    end
    
    % The following data may not used in the script, but will be used in
    % simulations. So, do not delete!
    [GmObj_Cell{i},GmDSS_Cell{i},ApparatusPara{i},ApparatusEqui{i},ApparatusDiscreDamping{i},OtherInputs{i},ApparatusStateStr{i},ApparatusInputStr{i},ApparatusOutputStr{i}] = ...
        SimplusGT.Toolbox.ApparatusModelCreate(ApparatusBus{i},ApparatusType{i},ApparatusPowerFlow{i},Para{i},Ts,ListBusNew);
    x_e{i} = ApparatusEqui{i}{1};
    u_e{i} = ApparatusEqui{i}{2};
end

% ### Get the appended model of all apparatuses
fprintf('Getting the appended descriptor state space model of all apparatuses...\n')
GmObj = SimplusGT.Toolbox.ApparatusModelLink(GmObj_Cell);

% ### Get the model of whole system
fprintf('Getting the descriptor state space model of whole system...\n')
[GsysObj,GsysDSS,Port_v,Port_i,BusPort_v,BusPort_i] = ...
    SimplusGT.Toolbox.ConnectGmZbus(GmObj,ZbusObj,N_Bus);

% ### Whole-system admittance model
YsysObj = SimplusGT.ObjTruncate(GsysObj,Port_i,Port_v);
YsysDSS = YsysObj.GetDSS(YsysObj);   

% ### Chech if the system is proper
fprintf('Checking if the whole system is proper:\n')
if isproper(GsysDSS)
    fprintf('Proper!\n');
    fprintf('Calculating the minimum realization of the system model for later use...\n')
    % GminSS = minreal(GsysDSS);
    GsysSS = SimplusGT.dss2ss(GsysDSS);
    % This "minreal" function only changes the element sequence of state
    % vectors, but does not change the element sequence of input and output
    % vectors.
    InverseOn = 0;
else
    error('Error: System is improper, which has more zeros than poles.')
end
if SimplusGT.is_dss(GsysSS)
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
fprintf('Whole system port model (system object form): GsysObj\n')
fprintf('Whole system port model (descriptor state space form): GsysDSS\n')
if Enable_PrintOutput
    [SysStateString,SysInputString,SysOutputString] = GsysObj.GetString(GsysObj);
    fprintf('Print ports of GsysDSS:\n')
    SimplusGT.Toolbox.PrintSysString(ApparatusBus,ApparatusType,GmObj_Cell,ZbusObj);
	fprintf('Print power flow result:\n')
    fprintf('The format below is "| bus | P | Q | V | angle | omega |". P and Q are in load convention.\n')
    ListPowerFlow
end

fprintf('Other models: \n')
fprintf('Whole system port model (state space form): GminSS\n')
fprintf('Whole system admittance model (system object form): YsysObj\n')
fprintf('Whole system admittance model (descriptor state space form): YsysDSS\n')

%%
% ==================================================
% Check Stability
% ==================================================

fprintf('\n')
fprintf('==================================\n')
fprintf('Check Stability\n')
fprintf('==================================\n')

fprintf('Calculatting pole/zero...\n')
% pole_sys_ = pole(GsysDSS)/2/pi;
[~,EigenValueSys] = eig(GsysSS.A);
EigenValueSys = diag(EigenValueSys);
EigenValueSys = EigenValueSys(find(real(EigenValueSys) ~= inf));
EigenValueSys = EigenValueSys/2/pi;
fprintf('Checking if the system is stable:\n')
if isempty(find(real(EigenValueSys)>1e-6, 1))
    fprintf('Stable!\n');
else
    fprintf('Warning: Unstable!\n')
end


    
