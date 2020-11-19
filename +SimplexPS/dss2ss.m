% This function converts a dss system to ss system

% Author(s): Yitong Li

function Gss = dss2ss(Gdss) 

A = Gdss.A;
B = Gdss.B;
C = Gdss.C;
D = Gdss.D;

Gss = ss(A,B,C,D);

end