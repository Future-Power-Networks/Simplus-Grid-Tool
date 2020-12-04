% This function achieves "dss_SwitchInOut" when the input is a customized
% matlab system object. 

% Author(s): Yitong Li, Yunjie Gu

%% Notes:
%
% The state string will be modified in this function.

%%
function obj_new = ObjSwitchInOut(obj,length_sw)

% Load data
[~,G] = obj.GetDSS(obj);
[StateStr,InputStr,OutputStr] = obj.GetString(obj);

% Create a new object
obj_new =SimplexPS.Class.ModelBase;

% Switch
G_new = SimplexPS.DssSwitchInOut(G,length_sw);
obj_new.SetDSS(obj_new,G_new);

% Get the string
StateStr_new = StateStr;
InputStr_new = InputStr;
OutputStr_new = OutputStr;
for i = 1:length_sw
    % Add xi_i into state string vector, which is caused by calling
    % "dss_SwitchInOut"
    StateStr_new = [StateStr_new,strcat('xi_',num2str(i))];
    
    % Switch the corresponding input and output string
    InputStr_new{i} = OutputStr{i};
    OutputStr_new{i} = InputStr{i};
end
obj_new.SetString(obj_new,StateStr_new,InputStr_new,OutputStr_new);

% Check the dimension
SimplexPS.ObjCheckDim(obj_new);

end