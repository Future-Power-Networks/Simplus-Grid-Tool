% Convert a state-space model to a symbolic transfer function

% Author(s): Yunjie Gu

%%
% this is very time consuming and should be used with care
% 2020-02-09 accelarated by diagnolizing A, much faster but still slow
% the accelarated version is acceptable and used to replace tf2sym to solve
% overflow problem of tf

%%
function Gsym = ss2sym(Gss,nu,ny)

    % Check if the system is in dss form or ss form.
    if SimplusGT.is_dss(Gss)
        error(['Error: The system is in dss form.']);
    end

    try 
        nu; %#ok<VUNUS>
    catch
        nu = 1:length(Gss.B(1,:));
    end
    
    try 
        ny; %#ok<VUNUS>
    catch
        ny = 1:length(Gss.C(:,1));
    end
    
    s = sym('s');
    
    [V,D] = eig(Gss.A);
    e = diag(D);
    h = 1./vpa(s*ones(length(e),1)-e);
    
    Gsym = zeros(length(ny),length(nu));
    Gsym = vpa(Gsym);
    for n = nu
        for m = ny
            b = V^(-1)*Gss.B(:,n);
            b = b.';
            c = Gss.C(m,:)*V;
            f = b.*c;
            Gsym(m,n) = vpa(f)*h + Gss.D(m,n);
        end
    end

    %I = eye(length(S.A));
    %G = S.C(ny,:)*(s*I-S.A)^(-1)*S.B(:,nu) + S.D(ny,nu);
    %H = diag(h);
    %G = S.C(ny,:)*V*H*V^(-1)*S.B(:,nu) + S.D(ny,nu);    
end