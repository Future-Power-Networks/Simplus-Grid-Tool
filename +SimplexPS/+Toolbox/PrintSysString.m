% This function prints the strings of the obtained state space system. 

% Author(s): Yitong Li

function PrintSysString(N_Device,DeviceType,DeviceStateStr,DeviceInputStr,DeviceOutputStr,ZbusStateStr)

%% Print state string
fprintf('\n')
fprintf('  Model state in order:\n');
CountState = 0;
% Device
IndexState{1} = 1;
for i = 1:N_Device
    fprintf(['    Device',num2str(i),':\n']);
    IndexState{i+1} = SimplexPS.PrintIndexCell(DeviceStateStr{i},6,IndexState{i});
    IndexState{i+1} = IndexState{i+1} + 1;
end
% Network
fprintf(['    Network line:\n']);
SimplexPS.PrintIndexCell(ZbusStateStr,6,IndexState{N_Device+1});

%% Print input string
fprintf('\n')
fprintf('  Model input in order:\n')
% Print device input string
IndexInput{1} = 1;
for i = 1:N_Device
    fprintf(['    Bus',num2str(i),':\n'])
    IndexInput{i+1} = SimplexPS.PrintIndexCell(DeviceInputStr{i},6,IndexInput{i});
    IndexInput{i+1} = IndexInput{i+1} + 1;
end

%% Print output string
fprintf('\n')
fprintf('  Model output in order:\n');
% Print device input string
IndexOutput{1} = 1;
for i = 1:N_Device
    fprintf(['    Bus',num2str(i),':\n']);
    IndexOutput{i+1} = SimplexPS.PrintIndexCell(DeviceOutputStr{i},6,IndexOutput{i});
    IndexOutput{i+1} = IndexOutput{i+1} + 1;
end

fprintf('\n')
end