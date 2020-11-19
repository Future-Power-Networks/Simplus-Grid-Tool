% This function gets the appended descriptor (implicit) state space system
% from two systems G1 and G2.

% Author(s): Yitong Li

%% Notes
%
% G1 and/or G2 can be empty system, static system.
%
% Example 2 can be seen as a special case of example 1, which can also be
% dealt with by this function, because the empty matrix also keeps
% dimensions (for example 0*2 empty matrix).

%% Example

% Example 1:
% G1
% E*dx1/dt = A1*x1 + B1*u1
% y1       = C1*x1 + D1*u1
% G2
% E*dx2/dt = A2*x2 + B2*u2
% y2       = C2*x2 + D2*u2
%
% -> Gss
% State equation
% [E1 0 ]*[dx1]/dt = [A1 0 ]*[x1] + [B1 0 ]*[u1]
% [0  E2] [dx2]      [0  A2] [x2]   [0  B2] [u2]
% Output equation
% [y1] = [C1 0 ]*[x1] + [D1 0 ]*[u1]
% [y2]   [0  C2] [x2]   [0  D2] [u2]

% Example 2:
% If G1 is a static system with A1=B1=C1=[]
%
% -> Gss
% Station equation
% E2*dx2/dt = [A2]*[x2] + [0 B2]*[u1]
%                                [u2]
% Output equation
% [y1] = [0 ]*[x2] + [D1 0 ]*[u1]
% [y2]   [C2]        [0  D2] [u2]

%%
function Gdss = DssAppend(G1,G2)

% Get the matrices
A1 = G1.A; B1 = G1.B; C1 = G1.C; D1 = G1.D; E1 = G1.E;
A2 = G2.A; B2 = G2.B; C2 = G2.C; D2 = G2.D; E2 = G2.E;

% Get the length of input, output, state
lx1 = length(A1); [ly1,lu1] = size(D1);
lx2 = length(A2); [ly2,lu2] = size(D2);

% Check if descriptor form
if ( isempty(E1) && (lx1~=0) ); E1 = eye(lx1); end
if ( isempty(E2) && (lx2~=0) ); E2 = eye(lx2); end

% Append
A = blkdiag(A1,A2);
B = blkdiag(B1,B2);
C = blkdiag(C1,C2);
E = blkdiag(E1,E2);
D = blkdiag(D1,D2);

% Get the new system    
Gdss = dss(A,B,C,D,E);

end