% This function constructs the summed descriptor (implicit) state space
% system G = G1 + G2.

% Author(s): Yitong Li

%% Notes

% Be very careful about the order (sequence) of elements in inputs and
% outputs vectors of two systems. For instance, if the inputs of G1 and
% G2 are [vd,vq] and [vq,vd] respectively, the calculated results would
% be WRONG.

% G1 and G2 can be static stata space model, which has A=B=C=E=[] and
% D~=[], as long as the inputs (and outputs) are in order.

% G1 and G2 CAN NOT be empty systems, because this means inputs (and/or
% outputs) are not in order.

%% Example

% ### Example 1: G1 and G2 have the same input and output
%
% Transfer function
% y = G*u = (G1 + G2)*u
% where G1 and G2 are two descriptor state space systems
%
% -> whole system state space:
% State equations
% [E1   ]*[dx1]/dt = [A1   ]*[x1] + [B1]*u
% [   E2] [dx2]      [   A2] [x2]   [B2]
% Output equations
% y = [C1 C2]*[x1] + (D1+D2)*u
%             [x2]

% ### Example 2: G1 and G2 have different input and output
%
% (Example 2 is not that genenral, it is better to use "append", "sumin",
% "sumout" to get a more general algorithm)
% 
% G1:
% with normal length of u and y
% State equation
% E1*dx1/dt = A1*x1 + [B1.a B1.b]*[u ]
%                                 [u']
% Output equation
% [y1 ] = [C1.a]*x1 + [D1.a D1.b]*[u ]
% [y1']   [C1.c]      [D1.c D1.d] [u']
% G2:
% State equation
% E2*dx2/dt = A2*x2 + B2*u
% Output equation
% y2 = C2*x2 + D2*u
%
% If [y ] = [y1+y2] -> whole system state space model:
%    [y']   [y1'  ]
% State equation
% [E1   ]*[dx1]/dt = [A1   ]*[x1] + [B1.a B1.b]*[u ]
% [   E2] [dx2]      [   A2] [x2]   [B2       ] [u']
% Output equation
% [y ]=[C1.a C2]*[x1] + [D1.a+D2 D1.b]*[u ]
% [y'] [C1.c   ] [x2]   [D1.c    D1.d] [u']

%% Function

function G = dss_Sum(G1,G2)

A1 = G1.A; B1 = G1.B; C1 = G1.C; D1 = G1.D; E1 = G1.E;
A2 = G2.A; B2 = G2.B; C2 = G2.C; D2 = G2.D; E2 = G2.E; 

% Get the length of state, input, output 
[lx1,lu1,ly1] = dss_GetDim(G1);
[lx2,lu2,ly2] = dss_GetDim(G2);

% Check for implicit state space
if isempty(E1); E1 = eye(lx1); end
if isempty(E2); E2 = eye(lx2); end

% Do the calculation
if ( (lu1==lu2) && (ly1==ly2) )
    A = blkdiag(A1,A2);
    B = [B1;B2];
    C = [C1,C2];
    D = D1+D2;
    E = blkdiag(E1,E2);
% elseif ( (lu1>lu2) && (ly1>ly2) )
%     % Assume:
%     % Example 2 is used
%     % [u;u'] with common input u 
%     % [y;y'] with common output y
%     % If input is [u';u] or output is [y';y], this function will get
%     % the WRONG results.
% 
%     % Cut the state space of system 1
%     B1_a = B1(:,1:lu2);
%     B1_b = B1(:,(lu2+1):end);
% 
%     C1_a = C1(1:ly2,:);
%     C1_c = C1((ly2+1):end,:);
% 
%     D1_a = D1(1:ly2,1:lu2);
%     D1_b = D1(1:ly2,(lu2+1):end);
%     D1_c = D1((ly2+1):end,1:lu2);
%     D1_d = D1((ly2+1):end,(lu2+1):end);
% 
%     % Calculate the zero matrix
%     [r_B2,~]   = size(B2);
%     [~,c_B1_b] = size(B1_b);
%     B0 = zeros(r_B2,c_B1_b);
% 
%     [~,c_C2]   = size(C2);
%     [r_C1_c,~] = size(C1_c);
%     C0 = zeros(r_C1_c,c_C2);
% 
%     % Form the new state space
%     A = blkdiag(A1,A2);
%     B = [B1_a, B1_b;
%          B2,   B0];
%     C = [C1_a, C2;
%          C1_c, C0];
%     D = [D1_a+D2, D1_b;
%          D1_c,    D1_d];
%     E = blkdiag(E1,E2);
else
    error(['Error']);
end

G = dss(A,B,C,D,E);
        
end