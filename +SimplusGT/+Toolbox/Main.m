% This is the main function for SimplusGT.

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
fprintf('Start to run Simplus Grid Tool\n')
fprintf('==================================\n')

%% 
% ==================================================
% Change current path of matlab
% ==================================================
PathStr = mfilename('fullpath');        % Get the path of main.m
[PathStr,~,~]  = fileparts(PathStr);    % Get the path of Toolbox namespace
[PathStr,~,~]  = fileparts(PathStr);    % Get the of root namespace
[PathStr,~,~]  = fileparts(PathStr);    % Get the path of toolbox
cd(PathStr);                            % Change the current address

%%
% ==================================================
% Load customized data
% ==================================================
fprintf('Loading data, please wait a second...\n')

% ### Re-arrange basic settings
Fs = InputData.Basic.Fs;
Ts = 1/Fs;               % (s), sampling period
Fbase = InputData.Basic.Fbase; % (Hz), base frequency
Sbase = InputData.Basic.Sbase; % (VA), base power
Vbase = InputData.Basic.Vbase; % (V), base voltage
Ibase = Sbase/Vbase;     % (A), base current
Zbase = Vbase/Ibase;     % (Ohm), base impedance
Ybase = 1/Zbase;         % (S), base admittance
Wbase = Fbase*2*pi;      % (rad/s), base angular frequency
Advance = InputData.Advance;
% Notes:
% The base values would be used in simulations, and should not be deleted
% here.

% Initialize figure index
FigN = 1000;

% ### Re-arrange the bus netlist
[ListBus,N_Bus] = SimplusGT.Toolbox.RearrangeListBusStruct(InputData);

% ### Re-arrange the line netlist
[ListLine,N_Branch.N_Bus_] = SimplusGT.Toolbox.RearrangeListLineStruct(InputData,ListBus);
DcAreaFlag = find(ListBus(:,12)==2);

% ### Re-arrange the apparatus netlist
NumApparatus = length(InputData.Apparatus);
for i = 1:NumApparatus
    ApparatusBus{i} = InputData.Apparatus(i).BusNo;
    ApparatusType{i} = InputData.Apparatus(i).Type;
    Para{i} = InputData.Apparatus(i).Para;
end
% The names of "ApparatusType" and "Para" can not be changed, because they
% will also be used in simulink model.

% Notes:
% No error checking if number of apparatuses is different from number of buses.

%%
% ==================================================
% Power flow analysis
% ==================================================

% ### Power flow analysis
fprintf('Do the power flow analysis...\n')
if ~isempty(DcAreaFlag)
    InputData.Advance.PowerFlowAlgorithm = 1;
    fprintf(['Warning: Because the system has dc area(s), the Gauss-Seidel power flow method is always used.\n']);
end
switch InputData.Advance.PowerFlowAlgorithm
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

% Move load flow (PLi and QLi) to bus admittance matrix
[ListBusNew,ListLineNew,PowerFlowNew] = SimplusGT.PowerFlow.Load2SelfBranch(ListBus,ListLine,PowerFlow);

% For print
ListPowerFlow = SimplusGT.PowerFlow.Rearrange(PowerFlow);
ListPowerFlowNew = SimplusGT.PowerFlow.Rearrange(PowerFlowNew);

fprintf('Print power flow result:\n')
fprintf('The format below is "| bus | P | Q | V | angle | omega |". P and Q are in load convention.\n')
ListPowerFlow

%%
% ==================================================
% Descriptor state space model
% ==================================================

EnableStateSpaceModel = 1;
if EnableStateSpaceModel

% ### Get the model of lines
fprintf('Get the descriptor state space model of network lines...\n')

[YbusObj,YbusDSS,~] = SimplusGT.Toolbox.YbusCalcDss(ListBusNew,ListLineNew,Wbase);
% if Enable_RemoveDuplicateStates == 1
%     YbusObj = SimplusGT.Toolbox.RemoveYbusDuplicateStates(YbusObj);
% end
ZbusObj = SimplusGT.ObjSwitchInOut(YbusObj,length(YbusDSS));

% ### Get the models of bus apparatuses
fprintf('Get the descriptor state space model of bus apparatuses...\n')
for i = 1:NumApparatus
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
fprintf('Get the appended descriptor state space model of all apparatuses...\n')
GmObj = SimplusGT.Toolbox.ApparatusModelLink(GmObj_Cell);

% ### Get the model of whole system
fprintf('Get the descriptor state space model of whole system...\n')
[GsysObj,GsysDSS,Port_v,Port_i,BusPort_v,BusPort_i] = ...
    SimplusGT.Toolbox.ConnectGmZbus(GmObj,ZbusObj,N_Bus);

% ### Whole-system admittance model
YsysObj = SimplusGT.ObjTruncate(GsysObj,Port_i,Port_v);
YsysDSS = YsysObj.GetDSS(YsysObj);   

% ### Chech if the system is proper
fprintf('Check if the whole system is proper:\n')
if isproper(GsysDSS)
    fprintf('Proper!\n');
    fprintf('Calculate the minimum realization of the system model for later use...\n')
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

% ### Print
fprintf('\n')
fprintf('Print State-Space Model: \n')
fprintf('Whole system port model (system object form): GsysObj\n')
fprintf('Whole system port model (descriptor state space form): GsysDSS\n')
if InputData.Advance.EnablePrintOutput
    [SysStateString,SysInputString,SysOutputString] = GsysObj.GetString(GsysObj);
    fprintf('Print ports of GsysDSS:\n')
    SimplusGT.Toolbox.PrintSysString(ApparatusBus,ApparatusType,GmObj_Cell,ZbusObj);
end

fprintf('Other models: \n')
fprintf('Whole system port model (state space form): GminSS\n')
fprintf('Whole system admittance model (system object form): YsysObj\n')
fprintf('Whole system admittance model (descriptor state space form): YsysDSS\n')

% ### Check stability
fprintf('\n')
fprintf('Calculate pole/zero...\n')
% pole_sys_ = pole(GsysDSS)/2/pi;
[~,EigenValueSys] = eig(GsysSS.A);
EigenValueSys = diag(EigenValueSys);
EigenValueSys = EigenValueSys(find(real(EigenValueSys) ~= inf));
EigenValueSys = EigenValueSys/2/pi;
fprintf('Check if the system is stable:\n')
if isempty(find(real(EigenValueSys)>1e-6, 1))
    fprintf('Stable!\n');
else
    fprintf('Warning: Unstable!\n')
end

% ### Plot fundamentals
fprintf('\n')
fprintf('Plot Fundamentals:\n')

% Plot pole/zero map
if InputData.Advance.EnablePlotPole
    fprintf('Plot pole map...\n')
    FigN = FigN+1;
    figure(FigN);
    subplot(1,2,1)
    scatter(real(EigenValueSys),imag(EigenValueSys),'x','LineWidth',1.5); hold on; grid on;
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
    title('Global pole map');
    
	subplot(1,2,2)
    scatter(real(EigenValueSys),imag(EigenValueSys),'x','LineWidth',1.5); hold on; grid on;
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
    title('Zoomed pole map');
    axis([-80,20,-150,150]);
    
    %SimplusGT.mtit('Pole Map');
else
    fprintf('Warning: The default plot of pole map is disabled.\n')
end

omega_p = logspace(-1,4,1e3)*2*pi;
omega_pn = [-flip(omega_p),omega_p];

% Plot admittance
if InputData.Advance.EnablePlotAdmittance
    fprintf('Plot admittance spectrum...')
  	FigN = FigN+1;
 	figure(FigN);
    CountLegend = 0;
    VecLegend = {};
    for k = 1:N_Bus
        [k1,k2] = SimplusGT.CellFind(ApparatusBus,k);
        % Plot the active bus admittance only
        if (0<=ApparatusType{k2} && ApparatusType{k2}<90) || ...
           (1000<=ApparatusType{k2} && ApparatusType{k2}<1090) || ...
           (2000<=ApparatusType{k2} && ApparatusType{k2}<2090)
           	Yss{k}  = minreal(GsysSS(BusPort_i{k},BusPort_v{k}));
            Ysym{k} = SimplusGT.ss2sym(Yss{k});
            SimplusGT.bode_c(Ysym{k}(1,1),1j*omega_p,'PhaseOn',0); 
            CountLegend = CountLegend + 1;
            VecLegend{CountLegend} = ['Bus',num2str(k)];
        end
    end
    legend(VecLegend);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (pu)');
    SimplusGT.mtit('Admittance Spectrum');
    
    save('Results\Yss','Yss');
else
    fprintf('Warning: The default plot of admittance spectrum is disabled.\n')
end

%
% ==================================================
% Modal Analysis
% ==================================================

fprintf('\n')
fprintf('==================================\n')
fprintf('Modal Analysis\n')
fprintf('==================================\n')
if (InputData.Advance.EnableParticipation == 1) && (isempty(DcAreaFlag))
    SimplusGT.Modal.ModalPreRun;
    SimplusGT.Modal.ModalAnalysis;
    fprintf('Generate GreyboxConfg.xlsx for user to config Greybox analysis.\n');    
else
    fprintf('Warning: The modal (participation) analysis is disabled or the power system has a dc area.\n');
end

else
    fprintf('Warning: The state space modeling is disabled.\n');
end

%%
% ==================================================
% Synchronization analysis
% ==================================================
fprintf('==================================\n')
fprintf('Synchronization Analysis\n')
fprintf('==================================\n')
EnableSynchronisationAnalysis = 0;
if EnableSynchronisationAnalysis
    SimplusGT.Synchron.MainSynchron();
else
    fprintf('Warning: The synchronisation analysis is disabled.\n')
end

%%
% ==================================================
% Create Simulink Model
% ==================================================
fprintf('\n')
fprintf('=================================\n')
fprintf('Simulink Model\n')
fprintf('=================================\n')

if N_Bus>=150
    InputData.Advance.EnableCreateSimulinkModel = 0;
    fprintf('Warning: The system has more than 150 buses;\n')
    fprintf('         The simulink model can not be created because of the limited size of GUI.\n')
    fprintf('         The static and dynamic analysis will not be influenced.\n')
end

if InputData.Advance.EnableCreateSimulinkModel == 1
    
    fprintf('Create the simulink model automatically, please wait a second...\n')

    % Set the simulink model name
    NameModel = 'mymodel_v1';

    % Close existing model with same name
    close_system(NameModel,0);
    
    % Create the simulink model
    SimplusGT.Simulink.MainSimulink(NameModel,ListBusNew,ListLineNew,ApparatusBus,ApparatusType,Advance);
    fprintf('Get the simulink model successfully! \n')
    fprintf('Please click the "run" button in the model to run it.\n')
    %fprintf('Warning: for later use of the simulink model, please "save as" a different name.\n')

else
    fprintf('Warning: The auto creation of simulink model is disabled.\n')
end

%%
fprintf('\n')
fprintf('==================================\n')
fprintf('End: run successfully.\n')
fprintf('==================================\n')