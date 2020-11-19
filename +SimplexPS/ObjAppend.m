% This function appends two matlab system objects

% Author(s): Yitong Li

%%
function Gobj = ObjAppend(Gobj1,Gobj2)

% Get the DSS models
[~,G1] = Gobj1.ReadDSS(Gobj1);
[~,G2] = Gobj2.ReadDSS(Gobj2);

% Append the models
G = SimplexPS.DssAppend(G1,G2);

% Get the strings
[StateStr1,InputStr1,OutputStr1] = Gobj1.ReadString(Gobj1);
[StateStr2,InputStr2,OutputStr2] = Gobj2.ReadString(Gobj2);

% Connect the strings
StateStr = [StateStr1,StateStr2];
InputStr = [InputStr1,InputStr2];
OutputStr = [OutputStr1,OutputStr2];

% Create a new object
Gobj = Class_Model_Base;
Gobj.LoadDSS(Gobj,G);
Gobj.WriteString(Gobj,StateStr,InputStr,OutputStr);

end