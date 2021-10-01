% This function achieves "dss_SwitchInOut" when the input is a customized
% matlab system object. 

% Author(s): Yitong Li, Yunjie Gu

%% Notes:
%
% The state string will be modified in this function.

%%
function ObjNew = ObjSwitchInOut(Obj,LengthSwtich)

% Load data
[~,G] = Obj.GetDSS(Obj);
[StateStr,InputStr,OutputStr] = Obj.GetString(Obj);

% Create a new object
ObjNew =SimplusGT.Class.ModelBase;

% Switch
G_New = SimplusGT.DssSwitchInOut(G,LengthSwtich);
ObjNew.SetDSS(ObjNew,G_New);

% Get the string
StateStrNew = StateStr;
InputStrNew = InputStr;
OutputStrNew = OutputStr;
for i = 1:LengthSwtich
    % Add xi_i into state string vector, which is caused by calling
    % "dss_SwitchInOut"
    StateStrNew = [StateStrNew,strcat('xi_',num2str(i))];
    
    % Switch the corresponding input and output string
    InputStrNew{i} = OutputStr{i};
    OutputStrNew{i} = InputStr{i};
end
ObjNew.SetString(ObjNew,StateStrNew,InputStrNew,OutputStrNew);

% Check the dimension
SimplusGT.ObjCheckDim(ObjNew);

end