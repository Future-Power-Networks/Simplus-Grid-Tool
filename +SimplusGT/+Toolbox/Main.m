% Main function for SimplusGT.

% Author(s): Yitong Li, Yunjie Gu
%
% Modified by Rob Oldaker, Yitong Li:
% # The function supports json user data input.
% # The function uses data by struct type.

%% Chech the input data type and convert it to struct
fprintf('\n')
fprintf('==================================\n')
fprintf('Check User Data File\n')
fprintf('==================================\n')
if UserDataType == 1
    UserData = [UserDataName, '.xlsx'];
    if isempty(which(UserData))
        UserData = [UserDataName,'.xlsm'];
        if isempty(which(UserData))
            error(['Error: The Excel format of "' UserDataName '" cannot be found.'])
        end
    end
    SimplusGT.Toolbox.Excel2Json(UserData);
elseif UserDataType == 0
else
    error('Error: Please check the setting of "UserDataType", which should be "1" or "0".')
end
UserData = [UserDataName, '.json'];
if isempty(which(UserData))
    error(['Error: The json format of "' UserDataName '" cannot be found.'])
end
fprintf([UserData,' is used for analysis.\n'])
which(UserData)
UserDataStruct = SimplusGT.JsonDecoder(UserData);

%%
fprintf('\n')
fprintf('==================================\n')
fprintf('Start: Run Simplus Grid Tool\n')
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
Fs = UserDataStruct.Basic.Fs;
Ts = 1/Fs;               % (s), sampling period
Fbase = UserDataStruct.Basic.Fbase; % (Hz), base frequency
Sbase = UserDataStruct.Basic.Sbase; % (VA), base power
Vbase = UserDataStruct.Basic.Vbase; % (V), base voltage
Ibase = Sbase/Vbase;     % (A), base current
Zbase = Vbase/Ibase;     % (Ohm), base impedance
Ybase = 1/Zbase;         % (S), base admittance
Wbase = Fbase*2*pi;      % (rad/s), base angular frequency
Advance = UserDataStruct.Advance;
% Notes:
% The base values would be used in simulations, and should not be deleted
% here.

% ### Re-arrange the bus netlist
[ListBus,NumBus] = SimplusGT.Toolbox.RearrangeListBusStruct(UserDataStruct);

% ### Re-arrange the line netlist
[ListLine,~,~] = SimplusGT.Toolbox.RearrangeListLineStruct(UserDataStruct,ListBus);
DcAreaFlag = find(ListBus(:,12)==2);

% ### Re-arrange the apparatus netlist
NumApparatus = length(UserDataStruct.Apparatus);
for i = 1:NumApparatus
    ApparatusBus{i} = UserDataStruct.Apparatus(i).BusNo;
    ApparatusType{i} = UserDataStruct.Apparatus(i).Type;
    Para{i} = UserDataStruct.Apparatus(i).Para;
end
clear('i');
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
    UserDataStruct.Advance.PowerFlowAlgorithm = 1;
    fprintf(['Warning: Because the system has dc area(s), the Gauss-Seidel power flow method is always used.\n']);
end
switch UserDataStruct.Advance.PowerFlowAlgorithm
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

[ObjYbusDss,YbusDss,~] = SimplusGT.Toolbox.YbusCalcDss(ListBusNew,ListLineNew,Wbase);
ObjZbusDss = SimplusGT.ObjSwitchInOut(ObjYbusDss,length(YbusDss));

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
    [ObjGmCell{i},GmDssCell{i},ApparatusPara{i},ApparatusEqui{i},ApparatusDiscreDamping{i},OtherInputs{i},ApparatusStateStr{i},ApparatusInputStr{i},ApparatusOutputStr{i}] = ...
        SimplusGT.Toolbox.ApparatusModelCreate(ApparatusBus{i},ApparatusType{i},ApparatusPowerFlow{i},Para{i},Ts,ListBusNew);
    x_e{i} = ApparatusEqui{i}{1};
    u_e{i} = ApparatusEqui{i}{2};
end
clear('i');

% ### Get the appended model of all apparatuses
fprintf('Get the appended descriptor state space model of all apparatuses...\n')
ObjGm = SimplusGT.Toolbox.ApparatusModelLink(ObjGmCell);

% ### Get the model of whole system
fprintf('Get the descriptor state space model of whole system...\n')
[ObjGsysDss,GsysDss,PortV,PortI,PortBusV,PortBusI] = ...
    SimplusGT.Toolbox.ConnectGmZbus(ObjGm,ObjZbusDss,NumBus);

% ### Whole-system admittance model
ObjYsysDss = SimplusGT.ObjTruncate(ObjGsysDss,PortI,PortV);
[~,YsysDss] = ObjYsysDss.GetDSS(ObjYsysDss); 
ObjYsysSs = SimplusGT.ObjDss2Ss(ObjYsysDss);
[~,YsysSs] = ObjYsysSs.GetSS(ObjYsysSs); 

% ### Chech if the system is proper
fprintf('Check if the whole system is proper:\n')
if isproper(GsysDss)
    fprintf('Proper!\n');
    fprintf('Calculate the minimum realization of the system model for later use...\n')
    ObjGsysSs = SimplusGT.ObjDss2Ss(ObjGsysDss);
    [~,GsysSs] = ObjGsysSs.GetSS(ObjGsysSs);
    % GsysMin = minreal(GsysSs);
else
    error('Error: GsysDss is improper, which has more zeros than poles.')
end

% ### Print
fprintf('\n')
fprintf('Print state space model: \n')
fprintf('Whole system port model (descriptor state space): GsysDss\n')
if UserDataStruct.Advance.EnablePrintOutput
    [GsysDssStateStr,GsysDssInStr,GsysDssOutStr] = ObjGsysDss.GetString(ObjGsysDss);
    [GsysSsStateStr,GsysSsInStr,GsysSsOutStr] = ObjGsysSs.GetString(ObjGsysSs);
    fprintf('Print ports of GsysDss:\n')
    SimplusGT.Toolbox.PrintSysString(ApparatusBus,ApparatusType,ObjGmCell,ObjZbusDss);
end

fprintf('Other models saved in workspace: \n')
fprintf('Whole system port model (state space): GsysSs\n')
fprintf('Whole system admittance model (descriptor state space): YsysDss\n')
fprintf('Whole system admittance model (state space): YsysSs\n')

% ### Check stability
fprintf('\n')
fprintf('Calculate pole/zero...\n')
[PhiMat,EigMat] = eig(GsysSs.A);
EigVec = diag(EigMat);
EigVec = EigVec(find(real(EigVec) ~= inf));
EigVecHz = EigVec/2/pi;
fprintf('Check if the system is stable:\n')
if isempty(find(real(EigVecHz)>1e-6, 1))
    fprintf('Stable!\n');
else
    fprintf('Warning: Unstable!\n')
end

% ### Plot fundamentals
fprintf('\n')
fprintf('Plot Fundamentals:\n')

% Plot pole/zero map
if UserDataStruct.Advance.EnablePlotPole
    fprintf('Plot pole map...\n')
    FigN = 100;
    PlotPoleMap(EigVecHz,FigN);
else
    fprintf('Warning: The default plot of pole map is disabled.\n')
end

% Plot admittance
if UserDataStruct.Advance.EnablePlotAdmittance
    fprintf('Plot admittance spectrum...\n')
  	FigN = 200;
    PlotAdmittanceSpectrum(NumBus,ApparatusBus,ApparatusType,GsysSs,PortBusI,PortBusV,FigN);
else
    fprintf('Warning: The default plot of admittance spectrum is disabled.\n')
end

%
% ==================================================
% Modal Analysis
% ==================================================

fprintf('\n')
fprintf('==================================\n')
fprintf('Modal Analysis: State Space \n')
fprintf('==================================\n')
if UserDataStruct.Advance.EnableParticipation == 1
    FigN = 300;
    ModeIndex = [1];
    SimplusGT.Modal.ModalAnalysisStateSpace(ObjGsysSs,ModeIndex,FigN);
else
    fprintf('Warning: This function is disabled.\n')
end

fprintf('\n')
fprintf('==================================\n')
fprintf('Modal Analysis: Transfer Function \n')
fprintf('==================================\n')
if 0
    SimplusGT.Modal.ModalPreRun;
    SimplusGT.Modal.ModalAnalysis;
    fprintf('Generate GreyboxConfg.xlsx for user to config Greybox analysis.\n');    
else
    fprintf('Warning: This function is disabled.\n');
end

else
    fprintf('Warning: The state space modeling is disabled.\n');
end

%%
% ==================================================
% Synchronization analysis
% ==================================================
fprintf('\n')
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

if NumBus>=150
    UserDataStruct.Advance.EnableCreateSimulinkModel = 0;
    fprintf('Warning: The system has more than 150 buses;\n')
    fprintf('         The simulink model can not be created because of the limited size of GUI.\n')
    fprintf('         The static and dynamic analysis will not be influenced.\n')
    fprintf('         This feature will be improved in the future version of SimplusGT.\n')
end

if UserDataStruct.Advance.EnableCreateSimulinkModel == 1
    
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
fprintf('End: Run Successfully.\n')
fprintf('==================================\n')

%%
% Author(s): Yitong Li

function PlotPoleMap(EigVecHz,FigN)
  
    figure(FigN);
    
    subplot(1,2,1)
    scatter(real(EigVecHz),imag(EigVecHz),'x','LineWidth',1.5); hold on; grid on;
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
    title('Global pole map');
    
	subplot(1,2,2)
    scatter(real(EigVecHz),imag(EigVecHz),'x','LineWidth',1.5); hold on; grid on;
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
    title('Zoomed pole map');
    axis([-80,20,-150,150]);
end

%%
% Author(s): Yitong Li

function PlotAdmittanceSpectrum(NumBus,ApparatusBus,ApparatusType,GsysSs,PortBusI,PortBusV,FigN)

    % Set frequency range 
    OmegaP = logspace(-1,4,500)*2*pi;
    OmegaPN = [-flip(OmegaP),OmegaP];
    CountLegend = 0;
    VecLegend = {};
    
    % Transform matrix from transfer function to complex vector
    T = [1,1i;
         1,-1i];
     
    for k = 1:NumBus
        [~,k2] = SimplusGT.CellFind(ApparatusBus,k);
        % Plot the active bus admittance only
        if (0<=ApparatusType{k2} && ApparatusType{k2}<90) || ...
           (1000<=ApparatusType{k2} && ApparatusType{k2}<1090) || ...
           (2000<=ApparatusType{k2} && ApparatusType{k2}<2090)
       
           	YcellSs{k}  = GsysSs(PortBusI{k},PortBusV{k});
            YcellSym{k} = SimplusGT.ss2sym(YcellSs{k});
            YcellSsCplx{k} = T*YcellSs{k}*T^(-1);
            YcellSymCplx{k} = SimplusGT.ss2sym(YcellSsCplx{k});
            
            figure(FigN);
            SimplusGT.bode_c(YcellSym{k}(1,1),1j*OmegaP,'PhaseOn',1); 
            figure(FigN+1);
            SimplusGT.bode_c(YcellSymCplx{k}(1,1),1j*OmegaPN,'PhaseOn',1); 
            
            CountLegend = CountLegend + 1;
            VecLegend{CountLegend} = ['Bus',num2str(k)];
        end
    end
 	figure(FigN)
  	SimplusGT.mtit('Transfer Function Matrix dq frame: Y_{dd}');
    legend(VecLegend);
  	figure(FigN+1)
  	SimplusGT.mtit('Complex Vector dq frame: Y_{dq+}');
    legend(VecLegend);
    
end


