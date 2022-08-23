% This function converts a dss system to ss system

% Author(s): Yitong Li 
% Modified by: Yunjie Gu

% Conditions to use this function: 
% 1: E must be diagonal

% format: input dss.(A,B,C,D,E), output ss.(A,B,C,D)
% if zeros > poles, output ss.(A,B,C,D,Bd,Dd) such that
% dx1/dt = Ax1 + Bu + Bd*du/dt
%      y = Cx1 + Du + Dd*du/dt
% where x1 is the states with non-zero E

% algorithms: 
% E dx/dt = Ax + Bu
% If E is diagonal and is not full rank, this equation can be re-written as
% E1 dx1/dt = A11*x1 + A12*x2 + B1*u (1)
%         0 = A21*x1 + A22*x2 + B2*u (2)
% if A22 is invertable, x2 can be solved from (2) directly
% if not, we find the maximum N so that N*A22 = 0 (N called null)
% and rank(N) + rank(A22) = length(A22)
% now we multiplies N to (2) to have
% N*A21*x1 + N*B2*u = 0              (3)
% N*A21*dx1/dt + N*B2*du/dt = 0      (4)
% N*A21*E1^(-1)[A11*x1 + A12*x2 + B1*u] + N*B2*du/dt = 0 
% define K = N*A21*E1^(-1)
% K * [A11*x1 + A12*x2 + B1*u] + N*B2*du/dt = 0    (5)
% combine (2) and (5)
% [ A21 ]*x1 + [ A22 ]*x2 + [ B2 ]*u + [  O ]*du/dt = 0
% [K*A11]      [K*A12]      [K*B1]     [N*B2] 
% or simply
%   A_21 *x1 +   A_22 *x2 +   B_2 *u +     F *du/dt = 0
% let W = left_inverse(A_22), that is,  W*A_22 = I (6)
% multiplies (5) by W
% x2 = -W*(A_21*x1 + B_2*u + F*du/dt)              (7)
% put (7) to (1)
% E1 dx1/dt = A11*x1 - A12*W*(A_21*x1 + B_2*u + F*du/dt) + B1*u  
% dx1/dt = E1^(-1)(A11 - A12*W*A_21) *x1   +   E1^(-1)(B1-A12*W*B_2) *u    +   -E1^(-1)A12*W*F *du/dt   (8)
%          ------------------------- new A     --------------------- new B     --------------- Bd
% y = C1*x1 + C2*x2 + Du
% y = C1*x1 + Du - C2*W*(A_21*x1 + B_2*u + F*du/dt)
% y = (C1-C2*W*A_21) *x1   +   (D-C2*W*B_2) *u   +   -C2*W*F *du/dt                                     (9)
%      ------------ new C       ---------- new D     ------- Dd

function [Gss,Index1] = dss2ss(Gdss) 

toler = 1e-14; % tolerance

A = Gdss.A;
B = Gdss.B;
C = Gdss.C;
D = Gdss.D;
E = Gdss.E;

if ~isdiag(E)
    error('Error: E matrix is not diagonal.');
end

DiagE = diag(E);
m = length(DiagE);
Index = zeros(m,1);
Point1 = 0;
Point2 = m+1;

% Find zero elements in E matrix
for n = 1:m
    if DiagE(n) < 1e-12             % Detect if E is zero with tolerance
        Point2 = Point2 - 1;
        Index(Point2) = n;          % Put zero element into E from E(end,end)
    else
        Point1 = Point1 + 1;
        Index(Point1) = n;          % Put non-zero element into E from E(1,1)
    end
end

Index1 = Index(1:Point1);    
Index2 = Index(Point2:m);

A11 = A(Index1,Index1);
A12 = A(Index1,Index2);
A21 = A(Index2,Index1);
A22 = A(Index2,Index2);

B1 = B(Index1,:);
B2 = B(Index2,:);

C1 = C(:,Index1);
C2 = C(:,Index2);

E1 = E(Index1,Index1);

% for i = 1:length(E1)
%     A11(i,:) = A11(i,:)/E1(i,i);
%     A12(i,:) = A12(i,:)/E1(i,i);
%     B1(i,:) = B1(i,:)/E1(i,i);
%     E1(i,i) = 1;
% end

if rank(A22,toler) == length(A22)
    InvA22 = A22^(-1);

    A_ = A11 - A12*InvA22*A21;
    B_ = B1  - A12*InvA22*B2;
    C_ = C1 - C2*InvA22*A21;
    D_ = D - C2*InvA22*B2;
    E_ = E1;

    InvE_ = E_^(-1);
    A_ = InvE_*A_;
    B_ = InvE_*B_;

    Gss = ss(A_,B_,C_,D_);
else
    N = null(A22.');
    N = N.';
    K = N*A21*E1^(-1);
    
    A_21 = [A21;K*A11];
    A_22 = [A22;K*A12];
    B_2  = [B2 ; K*B1];
    F    = [zeros(length(A22),length(B2(1,:)));N*B2];
    W    = SimplusGT.lft_inv(A_22);
    Ei   = E1^(-1);
    
    A_ = Ei*(A11 - A12*W*A_21);
    B_ = Ei*(B1-A12*W*B_2);
    Bd = -Ei*A12*W*F;
    
    C_ = C1-C2*W*A_21;
    D_ = D-C2*W*B_2;
    Dd = -C2*W*F;
    
    A_ = RemoveSmall(A_,toler);
    B_ = RemoveSmall(B_,toler);
    C_ = RemoveSmall(C_,toler);
    D_ = RemoveSmall(D_,toler);
    
    if (max(max(abs(Bd))) < toler) && (max(max(abs(Dd))) < toler)
        Gss = ss(A_,B_,C_,D_);
    else
        % note this is in struct not in ss form
        error(['Error: The system is improper.']);
        Gss.A = A_;
        Gss.B = B_;
        Gss.C = C_;
        Gss.D = D_;
        Gss.Bd = Bd;
        Gss.Dd = Dd;
    end
    
end
    
end

function X = RemoveSmall(Y,toler)
    X = Y;
    [N,M] = size(X);
    for m = 1:M
        for n = 1:N
            if(abs(X(n,m))<toler)
                X(n,m) = 0;
            end
        end 
    end
end