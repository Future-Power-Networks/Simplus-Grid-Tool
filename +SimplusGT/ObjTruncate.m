% This function truncates a state space system in object form.

% Author(s): Yitong Li

function [NewGobj] = ObjTruncate(Gobj,out,in)

% Get the model and strings of the old model
[~,Gdss] = Gobj.GetDSS(Gobj);
[StateStr,InStr,OutStr] = Gobj.GetString(Gobj);

% Update dss model
NewGdss = Gdss(out,in);

% Update strings
NewStateStr = StateStr;
for m = 1:length(in)
    NewInStr{m} = InStr{in(m)};
end
for n = 1:length(out)
    NewOutStr{n} = OutStr{out(n)};
end

% Create a new object model
NewGobj = SimplusGT.Class.ModelBase;
NewGobj.SetDSS(NewGobj,NewGdss);
NewGobj.SetString(NewGobj,NewStateStr,NewInStr,NewOutStr);
SimplusGT.ObjCheckDim(NewGobj);

end