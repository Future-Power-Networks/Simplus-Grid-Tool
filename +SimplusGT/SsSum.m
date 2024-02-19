% This function achieves the "sum" operation of two state space systems.

% Author(s): Yitong Li

function Gss = SsSum(Gss1, Gss2)

if is_dss(Gss1) || is_dss(Gss2)
    error(['Error: System 1 and/or 2 is not in ss form.']);
end

Gdss1 = ss2dss(Gss1);
Gdss2 = ss2dss(Gss2);

Gdss = dss_Sum(Gdss1,Gdss2);
Gss = dss2ss(Gdss);

end