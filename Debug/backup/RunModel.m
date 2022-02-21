%
% This script setsup a model prior to running simulink model
%
% Author(s): Rob Oldaker

modelFn = sprintf('Models/%s',modelName);
jsonFn = sprintf('Models/%s.json',modelName);
inputData = JsonDecoder(jsonFn);
%
% Setup model variables
%
SimplusGT.Toolbox.SetupModel();

%
% Run simulation
%
simOut = sim(modelFn,'SimulationMode','normal',...
            'SaveState','on','StateSaveName','xout',...
            'SaveOutput','on','OutputSaveName','yout',...
            'SaveFormat', 'Dataset');
%
% Save results to workspace variables
%
tout = simOut.tout;