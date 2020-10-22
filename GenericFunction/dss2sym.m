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
function Gsym = dss2sym(Gdss)

    A = Gdss.A;
    B = Gdss.B;
    C = Gdss.C;
    D = Gdss.D;
    E = Gdss.E;
    
    if isempty(E)
        error(['Error: The input system is not in dss form.']);
    end
    
    I = eye(length(E));
    s = sym('s');
    
    Gsym = C*inv(E*(s*I)-A)*B + D;
    
end