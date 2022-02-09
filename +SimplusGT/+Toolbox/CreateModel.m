%
% This script creates and saves a simulink model
%
% Author(s): Rob Oldaker
%

%
% Seup workspace variables
%
SimplusGT.Toolbox.SetupModel();
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

    
fprintf('Creating the simulink model automatically, please wait a second...\n')

% Set the simulink model name
Name_Model = 'mymodel_v1';

% Close existing model with same name
close_system(Name_Model,0);

% Create the simulink model
SimplusGT.Simulink.MainSimulink(Name_Model,ListBusNew,ListLineNew,ApparatusBus,ApparatusType,ListAdvance,PowerFlowNew);
[status,msg]=mkdir('Models'); % Returning "status" and "msg" stops it printing a warning
fn = strcat('Models/',modelName);
save_system(Name_Model,fn);
% output json as well
SaveAsJsonToFile(inputData,strcat(fn,'.json'));
close_system(fn)

fprintf('Simulink model %s successfully saved! \n', modelName);
%fprintf('Warning: for later use of the simulink model, please "save as" a different name.\n')

%%
fprintf('\n')
fprintf('==================================\n')
fprintf('End: run successfully.\n')
fprintf('==================================\n')

