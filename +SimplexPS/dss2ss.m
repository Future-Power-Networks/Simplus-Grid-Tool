% This function converts a dss system to ss system

% Author(s): Yitong Li 
% Modified by: Yunjie Gu

% Conditions to use this function: 
% 1: E must be diagonal
% 2: The system must have no more zeros than poles

function Gss = dss2ss(Gdss) 

A = Gdss.A;
B = Gdss.B;
C = Gdss.C;
D = Gdss.D;
E = Gdss.E;

if ~isdiag(E)
    error('E matrix is not diagonal');
end

DiagE = diag(E);
N = length(DiagE);
Index = zeros(N,1);
Point1 = 0;
Point2 = N+1;

for n = 1:N
    if DiagE(n) < 1e-12         %detect if E is zero with tolerance
        Point2 = Point2 - 1;
        Index(Point2) = n;
    else
        Point1 = Point1 + 1;
        Index(Point1) = n;        
    end
end

Index1 = Index(1:Point1);    
Index2 = Index(Point2:N);

A11 = A(Index1,Index1);
A12 = A(Index1,Index2);
A21 = A(Index2,Index1);
A22 = A(Index2,Index2);

B1 = B(Index1,:);
B2 = B(Index2,:);

C1 = C(:,Index1);
C2 = C(:,Index2);

E1 = E(Index1,Index1);

if rank(A22) < length(A22)
    error('The DSS system cannot be transformed to SS system. Try swap the input and output');
end
    
InvA22 = A22^(-1);

A_ = A11 - A12*InvA22*A21;
B_ = B1  - A12*InvA22*B2;
C_ = C1 - C2*InvA22*A21;
D_ = D;
E_ = E1;

InvE_ = E_^(-1);
A_ = InvE_*A_;
B_ = InvE_*B_;

Gss = ss(A_,B_,C_,D_);

end