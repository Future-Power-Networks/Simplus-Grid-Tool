% This function achieves "dss_SwitchInOut" when the input is a customized
% matlab system object. 

% Author(s): Yitong Li, Yunjie Gu

%% Notes:
%
% The state string will be modified in this function.

%%
function obj_new = obj_SwitchInOut(obj,length_sw)

% Load data
[~,G] = obj.ReadDSS(obj);
[StateStr,InputStr,OutputStr] = obj.ReadString(obj);

% Create a new object
obj_new = Class_Model_DSS;

% Switch
G_new = dss_SwitchInOut(G,length_sw);
obj_new.LoadModel(obj_new,G_new);

% Get the string
StateStr_new = StateStr;
InputStr_new = InputStr;
OutputStr_new = OutputStr;
for i = 1:length_sw
    % Add xi_i into state string vector, which is caused by calling
    % "dss_SwitchInOut"
    StateStr_new = [StateStr_new,strcat('xi_',num2str(i))];
end
obj_new.WriteString(obj_new,StateStr_new,InputStr_new,OutputStr_new);

% Check the dimension
obj_CheckDim(obj_new);

end