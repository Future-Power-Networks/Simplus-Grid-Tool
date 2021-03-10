% This function plots the Nyquist diagram of the symbolic system

% Author(s): Yitong Li

function [p1,p2] = nyquist_c(X,sbd,varargin)

% Set default values
Half = 0;
LineWidth = 1;
Color1 = 'k';
Color2 = 'k';

% Update values
for i = 1:length(varargin)
    if strcmpi(varargin{i},'Half')
        Half = varargin{i+1};
    elseif strcmpi(varargin{i},'LineWidth')
        LineWidth = varargin{i+1};
    elseif strcmpi(varargin{i},'Color1')
        Color1 = varargin{i+1};
 	elseif strcmpi(varargin{i},'Color2')
        Color2 = varargin{i+1};
    elseif strcmpi(varargin{i},'Color')
        Color1 = varargin{i+1};
        Color2 = varargin{i+1};
    end
end

funcX = matlabFunction(X);

Xw = zeros(1,length(sbd));

for n = 1:length(sbd)
    try 
        Xw(n) = funcX(sbd(n));
    catch
        Xw(n) = funcX();
    end
end

W = length(Xw);

if (imag(sbd(1))*imag(sbd(W)))<0
wbdn = sbd(1:W/2);
wbdp = sbd(W/2+1:W);

Xwn  = Xw(1:W/2);
Xwp  = Xw(W/2+1:W);
Xwn = flip(Xwn);

Arg_wn = angle(Xwn);
Arg_wp = angle(Xwp);

for n = 1:W/2
    Xwn_x(n) = real(Xwn(n));
    Xwn_y(n) = imag(Xwn(n));
    Xwp_x(n) = real(Xwp(n));
    Xwp_y(n) = imag(Xwp(n));
end

p1 = plot(Xwn_x,Xwn_y,'LineWidth',LineWidth,'Color',Color1); grid on; hold on;

if Half == 0
    p2 = plot(Xwp_x,Xwp_y,'LineWidth',LineWidth,'Color',Color2); grid on; hold on;
end

else
    error(['Error: Please add both positive and negative frequency']);
end

end