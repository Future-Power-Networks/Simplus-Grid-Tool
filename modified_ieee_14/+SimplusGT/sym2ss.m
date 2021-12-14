% This function converts a symbolic system to state space system

% Author(s): Yitong Li

% Notes:
%
% The state space system can not be in descriptor form.

function Y = sym2ss(X)

    [xn,xd] = numden(X);
    xnp = sym2poly(xn);
    xdp = sym2poly(xd);
    [A,B,C,D] = tf2ss(xnp,xdp);
    Y = ss(A,B,C,D);

end