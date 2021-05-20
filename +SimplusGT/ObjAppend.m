% This function appends two matlab system objects

% Author(s): Yitong Li

%%
function Gobj = ObjAppend(Gobj1,Gobj2)

% Get the DSS models
[~,G1] = Gobj1.GetDSS(Gobj1);
[~,G2] = Gobj2.GetDSS(Gobj2);

% Append the models
G = SimplusGT.DssAppend(G1,G2);

% Get the strings
[StateStr1,InputStr1,OutputStr1] = Gobj1.GetString(Gobj1);
[StateStr2,InputStr2,OutputStr2] = Gobj2.GetString(Gobj2);

% Connect the strings
StateStr = [StateStr1,StateStr2];
InputStr = [InputStr1,InputStr2];
OutputStr = [OutputStr1,OutputStr2];

% Create a new object
Gobj = SimplusGT.Class.ModelBase;
Gobj.SetDSS(Gobj,G);
Gobj.SetString(Gobj,StateStr,InputStr,OutputStr);

end