% convert a descriptor (implicit) state-space model to a symbolic transfer function

% Author(s): Yitong Li

%%
% Notes:
%
% This function is very time consuming, and should be used very carefully.
%
% Original system format
% Edx/dt = Ax + Bu
% y      = Cx + Du

%%
function Gsym = dss2sym(Gss)

    A = Gss.A;
    B = Gss.B;
    C = Gss.C;
    D = Gss.D;
    E = Gss.E;
    I = eye(length(E));
    s = sym('s');
    
    Gsym = C*inv(E*(s*I)-A)*B + D;  
end