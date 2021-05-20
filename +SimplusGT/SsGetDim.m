% This function gets the dimension of a state space system.

% Author(s): Yitong Li

function [lx,lu,ly] = SsGetDim(Gss)

if SimplusGT.is_dss(Gss)
    error(['Error: The system is not in ss form.']);
end

Gdss = SimplusGT.ss2dss(Gss);
[lx,lu,ly] = SimplusGT.DssGetDim(Gdss);

end