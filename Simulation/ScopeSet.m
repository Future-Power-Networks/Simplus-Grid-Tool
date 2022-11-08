% This function sets the scope data in the simulink model. So that the
% scope data can be saved into workspace.

% Author(s): Yitong Li

% ModelName = 'Sim68BusSg';
ModelName = 'Sim68BusCap';
% ModelName = 'Sim68BusUnbalancedFault';
% ModelName = 'Sim68BusGFL';

Fbus = [19,22,30,31,32,34,35,37,38,43,54,57,58,62,63,65,66];

% NumBus = 68;
NumBus = 16;
% NumBus = 26;
for i = 1:NumBus
    if isempty(find(Fbus == i,1))
        myConfiguration = get_param([ModelName '/D-Scope' num2str(i)],'ScopeConfiguration');
        myConfiguration.DataLogging = true;
        myConfiguration.DataLoggingDecimateData = true;
        % myConfiguration.DataLoggingDecimation = '30';       % For high inertia test
        myConfiguration.DataLoggingDecimation = '60';    % For low and medieum inertia test
        myConfiguration.DataLoggingVariableName = ['Data_App' num2str(i) ];
        myConfiguration.DataLoggingSaveFormat = 'Dataset';
    end
end

