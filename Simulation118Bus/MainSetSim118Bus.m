% Author(s): Yitong Li

clear all
clc
close all

ModelName = 'Sim118Bus';

Fbus = [2,3,5,7,9,11,13,14,16,17,20:23,28:30,33,35,37:39,41,43:45,47,48,50:53,57,58,60,63,64,67,68,71,75,78,79,81:84,86,88,93:98,101,102,106,108,109,114,115,117,118];

NumBus = 118;
% NumBus = 16;
for i = 1:NumBus
    if isempty(find(Fbus == i,1))
        myConfiguration = get_param([ModelName '/D-Scope' num2str(i)],'ScopeConfiguration');
        myConfiguration.DataLogging = true;
        myConfiguration.DataLoggingDecimateData = true;
        myConfiguration.DataLoggingDecimation = '30';
        myConfiguration.DataLoggingVariableName = ['App' num2str(i) ];
        myConfiguration.DataLoggingSaveFormat = 'Dataset';
    end
end
