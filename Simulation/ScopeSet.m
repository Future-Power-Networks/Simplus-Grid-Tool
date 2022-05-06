% Set scope

% ModelName = 'K_Sim';
% ModelName = 'K_Sim_NoCap';
% ModelName = 'Sim68BusSg';
ModelName = 'Sim68BusCap';

Fbus = [19,22,30,31,32,34,35,37,38,43,54,57,58,62,63,65,66];

% NumBus = 68;
NumBus = 16;
for i = 1:NumBus
    if isempty(find(Fbus == i,1))
        myConfiguration = get_param([ModelName '/D-Scope' num2str(i)],'ScopeConfiguration');
        myConfiguration.DataLogging = true;
        myConfiguration.DataLoggingDecimateData = true;
        myConfiguration.DataLoggingDecimation = '30';
        % myConfiguration.DataLoggingDecimation = '120';
        myConfiguration.DataLoggingVariableName = ['Data_App' num2str(i) ];
        myConfiguration.DataLoggingSaveFormat = 'Dataset';
    end
end

