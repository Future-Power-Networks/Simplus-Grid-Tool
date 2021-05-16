% This function plots the equivalent torque coefficient in complex plane.

% Author(s): Yitong Li

function plot_K(K,w,varargin)

[LineWidth,~] = LoadVar(1.5,'LineWidth',varargin);
[Color,~] = LoadVar([],'Color',varargin);

funcK = matlabFunction(K);

Kw = funcK(w);

if isempty(Color)
    scatter(real(Kw),imag(Kw),'.','LineWidth',LineWidth);
else
    scatter(real(Kw),imag(Kw),'.','LineWidth',LineWidth,'Color',Color);
end

end