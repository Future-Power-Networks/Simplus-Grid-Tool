%% case 1: 
%   ---[L1]-----[L3]-----[L5]--- +
%   |        |        |    <-- y
%   |       [C2]     [C4]        u
%   |        |        |
%   ---------------------------- -
%
%   L1=L5=1, C2=C4=L3=0

A = [ 0  1  0  0  0; ...
     -1  0  1  0  0; ...
      0 -1  0  1  0; ...
      0  0 -1  0  1; ...
      0  0  0 -1  0];
  
E = diag([1 0 0 0 1]);

C = [0 0 0 0 1];
B = C.';
D = 0;

Gdss.A = A;
Gdss.B = B;
Gdss.C = C;
Gdss.D = D;
Gdss.E = E;

Gss = SimplusGT.dss2ss(Gdss)

%% case 2:
Gdss = SimplusGT.DssSwitchInOut(Gdss,1);
Gss = SimplusGT.dss2ss(Gdss)
