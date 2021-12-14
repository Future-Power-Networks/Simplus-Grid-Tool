% This function converts a ss system to dss system

% Author(s): Yitong Li

function Gdss = ss2dss(Gss)

A = Gss.A;
B = Gss.B;
C = Gss.C;
D = Gss.D;
E = Gss.E;

if SimplusGT.is_dss(Gss)
    error(['Error: The original system is in dss system form.']);
end

E = eye(length(A));

Gdss = dss(A,B,C,D,E);

end