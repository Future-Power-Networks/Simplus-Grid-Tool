% This function gets the dimension of a state space system.

% Author(s): Yitong Li

function [lx,lu,ly] = SsGetDim(Gss)

if SimplexPS.is_dss(Gss)
    error(['Error: The system is not in ss form.']);
end

Gdss = SimplexPS.ss2dss(Gss);
[lx,lu,ly] = SimplexPS.DssGetDim(Gdss);

end