% This function gets the dimension of a state space system.

% Author(s): Yitong Li

function [lx,lu,ly] = ss_GetDim(Gss)

if is_dss(Gss)
    error(['Error: The system is not in ss form.']);
end

Gdss = ss2dss(Gss);
[lx,lu,ly] = dss_GetDim(Gdss);

end