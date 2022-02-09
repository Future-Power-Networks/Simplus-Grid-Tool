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
basicData = inputData.basic;
Fs = basicData.Fs;
Ts = 1/Fs;               % (s), sampling period
Fbase = basicData.Fbase; % (Hz), base frequency
Sbase = basicData.Sbase; % (VA), base power
Vbase = basicData.Vbase; % (V), base voltage
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

%%
% ==================================================
% Create Simulink Model
% ==================================================
fprintf('\n')
fprintf('=================================\n')
fprintf('Simulink Model\n')
fprintf('=================================\n')

if N_Bus>=150
    Enable_CreateSimulinkModel = 0;
    fprintf('Warning: The system has more than 150 buses;\n')
    fprintf('         The simulink model can not be created because of the limited size of GUI.\n')
    fprintf('         The static and dynamic analysis will not be influenced.\n')
end

if Enable_CreateSimulinkModel == 1
    
    fprintf('Creating the simulink model automatically, please wait a second...\n')

    % Set the simulink model name
    Name_Model = 'mymodel_v1';

    % Close existing model with same name
    close_system(Name_Model,0);
    
    % Create the simulink model
    SimplusGT.Simulink.MainSimulink(Name_Model,ListBusNew,ListLineNew,ApparatusBus,ApparatusType,ListAdvance,PowerFlowNew);
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
if Enable_PlotAdmittance
    fprintf('Plotting admittance spectrum...')
  	figure_n = figure_n+1;
 	figure(figure_n);
    CountLegend = 0;
    VecLegend = {};
    for k = 1:N_Bus
        [k1,k2] = SimplusGT.CellFind(ApparatusBus,k);
        % Plot the active bus admittance only
        if (0<=ApparatusType{k2} && ApparatusType{k2}<90) || ...
           (1000<=ApparatusType{k2} && ApparatusType{k2}<1090) || ...
           (2000<=ApparatusType{k2} && ApparatusType{k2}<2090)
           	Yss{k}  = GsysSS(BusPort_i{k},BusPort_v{k});
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
else
    fprintf('Warning: The default plot of admittance spectrum is disabled.\n')
end

%%
% ==================================================
% Modal Analysis
% ==================================================

fprintf('\n')
fprintf('==================================\n')
fprintf('Modal Analysis\n')
fprintf('==================================\n')
if (Enable_Participation == 1) && (isempty(DcAreaFlag))
    SimplusGT.Modal.ModalPreRun;
    SimplusGT.Modal.ModalAnalysis;
    fprintf('Generating GreyboxConfg.xlsx for user to config Greybox analysis.\n');    
else
    fprintf('Warning: The modal (participation) analysis is disabled or the power system has a dc area.');
end


%%
fprintf('\n')
fprintf('==================================\n')
fprintf('End: run successfully.\n')
fprintf('==================================\n')
%%
% Remained questions:
% - Advanced power flow, including droop bus, etc
% - The power flow calculation assumes the frequency is Wbase
% - Initialization of network lines, such as the current of line inductor
% and the voltage of line capacitor.
% - The calculation of discrete damping resistor should be double checked,
% especially for the interlink ac-dc converter which has hybrid ac-dc
% electrical ports.