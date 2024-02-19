% This function constructs a large state space model from a small state
% space model

% Author(s): Yitong Li

%% Example

% Example 1:
% Transfer function
% [y1] = [G11 G12]*[u1]
% [y2]   [G21 G22] [u2]
% where G11, G12, G21, G22 are four state space systems.
%
% -> whole system state space model:
% State equations
% [E11            ]*[dx11]/dt = [A11            ]*[x11] + [B11    ]*[u1]
% [    E12        ] [dx12]      [    A12        ] [x12]   [    B12] [u2]
% [        E13    ] [dx21]      [        A13    ] [x13]   [B21    ]
% [            E14] [dx22]      [            A14] [x14]   [    B22]
% Output equations
% [y1] = [C11 C12        ]*[x11] + [D11 D12]*[u1]
% [y2]   [        C21 C22] [x12]   [D21 D22] [u2]
%                          [x21]
%                          [x22]

%% Function

function G = DssArrange(Gcell)

[rmax,cmax] = size(Gcell);       % r -> row, c -> column

A = []; B = []; E = [];     % Reset
C = []; D = [];
for r = 1:rmax
    A1 = []; B1 = []; E1 = [];     % Reset
    C1 = []; D1 = [];

    % Get the state space model of one row
    for c = 1:cmax
        An = Gcell{r,c}.A;
        Bn = Gcell{r,c}.B;
        Cn = Gcell{r,c}.C;
        Dn = Gcell{r,c}.D;
        En = Gcell{r,c}.E;
        
        % Check if in desctiptor form
        if (isempty(En) && (~isempty(An)))
            En = eye(length(An));   % This will convert the ss form to dss form.
            % error(['Error: the system is not in descritpor-state-space form']);
        end
        
        % Re-arrange two cells horizontally
        A1 = blkdiag(A1,An);
        B1 = blkdiag(B1,Bn);
        C1 = [C1,Cn];
        D1 = [D1,Dn];
        E1 = blkdiag(E1,En);
    end

    % Re-arrange two cells virtically
    A = blkdiag(A,A1);
    B = [B;
         B1];
    C = blkdiag(C,C1);
    D = [D;
         D1];
    E = blkdiag(E,E1);
end

G = dss(A,B,C,D,E);

end