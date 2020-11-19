% This function switches the input and output of a descriptor (implicit)
% state space system by extending the state vector.

% Author(s): Yitong Li, Yunjie Gu

%%
% Notes:
%
% The key of this function is adding a new state xi, and letting the
% switched input u = xi.
%
% The system CAN be static system, but CAN NOT be empty system. In
% addition, the system has to have both input and output.
%
% The dimenstion of the switches input and output should be same. The heads
% of input and output are switched.
%
% This functon can deal with the appearance of not only the derivative but
% also the infinite gain.

%%
% Example 1:
% Original system:
% State equation
% E*dx/dt = A*x + [B1 B2]*[u1]
%                         [u2]
% Output equation
% [y1] = [C1]*x + [D11 D12]*[u1]
% [y2]   [C2]     [D21 D22] [u2]
%
% => switch u1 and y1 =>
% New system:
% State equation
% [E 0]*[dx ]/dt = [A   B1  ]*[x ] + [0 B2  ]*[y1]
% [0 0] [dxi]      [-C1 -D11] [xi]   [I -D12] [u2]
% Output equation
% [u1] = [0  I  ]*[x ] + [0 0  ]*[y1]
% [y2]   [C2 D21] [xi]   [0 D22] [u2]

% Example 2:
% If system is static, i.e., A, E, B, C are empty matrices
% Original system:
% State equation
% [] = []*[] + []*[u1]
%                 [u2]
% Output equation
% [y1] = [D11 D12]*[u1]
% [y2]   [D21 D22] [u2]
%
% => switch u1 and y1 =>
% New system:
% State equation
% [0]*[dxi]/dt = [-D11]*[xi] + [I -D12]*[y1]
%                                       [u2]
% Output equation
% [u1] = [I  ]*[xi] + [0 0  ]*[y1]
% [y2]   [D21]        [0 D22] [u2]

%%
function Gnew = DssSwitchInOut(G,lsw)

% Get the state space matrices
A = G.A;
B = G.B;
C = G.C;
D = G.D;
E = G.E;

% Get the length of input, output, state
lx = length(A);
[ly,lu] = size(D);
if ((lu==0) || (ly==0))
    error(['Errror: system has no input and/or output']);
elseif ((lsw>lu) || (lsw>ly))   % Check the switch range 
    error(['Error: switch range exceeds the length of input or output']); 
end

% Check if descriptor state space
if (isempty(E) && (lx~=0)); E = eye(lx); end 

% Cut B
B1 = B(:,1:lsw);
B2 = B(:,(lsw+1):end);

% Cut C
C1 = C(1:lsw,:);
C2 = C((lsw+1:end),:);

% Cut D
D11 = D(1:lsw,1:lsw);
D12 = D(1:lsw,(lsw+1):end);
D21 = D((lsw+1):end,1:lsw);
D22 = D((lsw+1):end,(lsw+1):end);

% For static system, A=B=C=[], B1=B2=[], C1=C2=[].
% For lsw=lu, B2=D12=D22=[].
% For lsw=ly, C2=D21=D22=[].

% Get the new matrices
Anew = [A,B1;
        -C1,-D11];
Bnew = [zeros(lx,lsw),B2;
        eye(lsw),-D12];
Cnew = [zeros(lsw,lx),eye(lsw);
        C2,D21];
Dnew = [zeros(lsw,lsw),zeros(lsw,lu-lsw);
        zeros(ly-lsw,lsw),D22]; 
Enew = blkdiag(E,zeros(lsw,lsw));

Gnew = dss(Anew,Bnew,Cnew,Dnew,Enew);

end