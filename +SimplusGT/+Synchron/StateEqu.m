% Simplified swing equations of a multi-apparatus power system
%
% Author(s): Yitong Li

function dxdt = StateEqu(x,GAMMA,gamma,Hmat,Dmat,Wref,Wbase)

Napp = length(Hmat);
H = diag(Hmat);
D = diag(Dmat);

% The state vector is represented as
% x = [theta1; theta2; ...; thetaN; omega1; omega2; ... omegaN]

% dtheta/dt = omega
for m = 1:Napp
    dxdt(m,1) = x(Napp+m);
end

% Swing equation
for m = 1:Napp
    K(m) = 0;
    for n = 1:Napp
        K(m) = K(m) + GAMMA(m,n)*sin(x(m) - x(n) + gamma(m,n));
    end
    % dxdt(m+Napp,1) = (Wref(m) - D(m)*(x(m+Napp)-Wbase) - K(m))/H(m);
    dxdt(m+Napp,1) = (Wref(m) - D(m)*(x(m+Napp)-Wbase) + K(m))/H(m);
    % dxdt(m+Napp,1) = (Wref(m) - D(m)*(x(m+Napp)) - K(m))/H(m);
end

end