% This function prints the strings of the obtained state space system. 

% Author(s): Yitong Li

%% Notes
%
% This function prints the string based on the order of apparatuses, because
% GmObj is obtained by appending the apparatus models.

%%

function PrintSysString(ApparatusBus,ApparatusType,ApparatusCell,ZbusObj)

%% Get strings
for n = 1:length(ApparatusCell)
    [ApparatusStateStr{n},ApparatusInStr{n},ApparatusOutStr{n}] = ApparatusCell{n}.GetString(ApparatusCell{n});
end

[ZbusStateStr,~,~] = ZbusObj.GetString(ZbusObj);

%%
N_Apparatus = length(ApparatusBus);

%% Print state string
fprintf('\n')
fprintf('  Model state in order:\n');
CountState = 0;
% Apparatus
IndexState{1} = 1;
for i = 1:N_Apparatus
    if ApparatusType{i} ~= 100 && ApparatusType{i} ~= 1100
        % Print only when the apparatus is not a floating bus.
        fprintf(['    Apparatus ',num2str(ApparatusBus{i}'),':\n']);
        IndexState{i+1} = SimplusGT.PrintIndexCell(ApparatusStateStr{i},6,IndexState{i});
        IndexState{i+1} = IndexState{i+1} + 1;
    else
        IndexState{i+1} = SimplusGT.PrintIndexCell(ApparatusStateStr{i},6,IndexState{i});
    end
end
% Network
fprintf(['    Network line:\n']);
SimplusGT.PrintIndexCell(ZbusStateStr,6,IndexState{N_Apparatus+1});

%% Print input string
fprintf('\n')
fprintf('  Model input in order:\n')
% Print apparatus input string
IndexInput{1} = 1;
for i = 1:N_Apparatus
    fprintf(['    Bus ',num2str(ApparatusBus{i}'),':\n'])
    IndexInput{i+1} = SimplusGT.PrintIndexCell(ApparatusInStr{i},6,IndexInput{i});
    IndexInput{i+1} = IndexInput{i+1} + 1;
end

%% Print output string
fprintf('\n')
fprintf('  Model output in order:\n');
% Print apparatus input string
IndexOutput{1} = 1;
for i = 1:N_Apparatus
    fprintf(['    Bus ',num2str(ApparatusBus{i}'),':\n']);
    IndexOutput{i+1} = SimplusGT.PrintIndexCell(ApparatusOutStr{i},6,IndexOutput{i});
    IndexOutput{i+1} = IndexOutput{i+1} + 1;
end

fprintf('\n')
end