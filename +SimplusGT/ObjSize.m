% Author(s): Yitong Li

function [LthState,LthIn,LthOut] = ObjSize(ObjG)

[StateStr,InStr,OutStr] = ObjG.GetString(ObjG);

LthState = length(StateStr);
LthIn = length(InStr);
LthOut = length(OutStr);

end