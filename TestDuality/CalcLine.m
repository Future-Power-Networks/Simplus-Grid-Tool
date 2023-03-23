% This function calculates the X or R from Z and X/R ratio.

% Author(s): Yitong Li

% Notes:

% Z can be an array of parameters.

function [X,R] = CalcLine(Z,XRratio)

for i = 1:length(Z)
    R(i) = sqrt(Z(i)^2/(1+XRratio^2));
    X(i) = R(i)*XRratio;
end

end