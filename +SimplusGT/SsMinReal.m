% Get the minimum realization of a state space model by Kalman
% decomposition.

% Author(s): Yitong Li

function GminKalman = SsMinReal(G)

    if (SimplusGT.is_dss(G))
        error('Error: The input system cannot be a descriptor state space model.')
    end
    A = G.A;
    B = G.B;
    C = G.C;
    D = G.D;
    
    [~,U] = minreal(G);
    A_ = U*A*(U');
    B_ = (U)*B;
    C_ = C*(U');
    D_ = D;
    
    Tor = 1e-9;
    A_ = TorCheck(A_,Tor);
    B_ = TorCheck(B_,Tor);
    C_ = TorCheck(C_,Tor);
    D_ = TorCheck(D_,Tor);
    
    GminKalman = ss(A_,B_,C_,D_);
    
end

function A = TorCheck(A,tor)
    [m,n] = size(A);
    for i = 1:m
        for j = 1:n
            if (A(i,j) <= tor)
                A(i,j) = 0;
            end
        end
    end
end