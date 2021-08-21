% convert a state-space model to a symbolic transfer function
% this is very time consuming and should be used with care
% 2020-02-09 accelarated by diagnolizing A, much faster but still slow
% the accelarated version is acceptable and used to replace tf2sym to solve
% overflow problem of tf

function G = ss2sym(S,nu,ny)

    try 
        nu; %#ok<VUNUS>
    catch
        nu = 1:length(S.B(1,:));
    end
    
    try 
        ny; %#ok<VUNUS>
    catch
        ny = 1:length(S.C(:,1));
    end
    
    s = sym('s');
    
    [V,D] = eig(S.A);
    e = diag(D);
    h = 1./vpa(s*ones(length(e),1)-e);
    
    G = zeros(length(ny),length(nu));
    G = vpa(G);
    for n = nu
        for m = ny
            b = V^(-1)*S.B(:,n);
            b = b.';
            c = S.C(m,:)*V;
            f = b.*c;
            G(m,n) = vpa(f)*h + S.D(m,n);
        end
    end

    %I = eye(length(S.A));
    %G = S.C(ny,:)*(s*I-S.A)^(-1)*S.B(:,nu) + S.D(ny,nu);
    %H = diag(h);
    %G = S.C(ny,:)*V*H*V^(-1)*S.B(:,nu) + S.D(ny,nu);    
end