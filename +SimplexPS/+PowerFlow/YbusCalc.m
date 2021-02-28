% This function calculates the "steady" nodal admittance matrix.

% Author(s): Yunjie Gu, Yitong Li

%% Notes:
%
% Format of branch:
%             ---             ---C---
% FromBus ---|a:1|---R---L---|       |--- ToBus
%             ---             ---G---
% where (a:1) is the turn ratio of transformer 
%
% Format of input:
% netlist:   |  From |  To   |   R   |  wL   |  wC   |   G   |
%            |  Bus  |  Bus  |       |       |       |       |
% linedata = [  1       2        0       X       0      inf;
%               1       1        0       0       B       G;
%               2       2        0       0       0       0];
%
% The obtained nodal admittance matrix is "steady", which means the
% obtained matrix does NOT contain Laplace operator "s".
%
% The matrix is in phasor frame.
%
% This function supports both ac and dc grid calculation. For dc grids,
% only resistance is considerred, because jX=0 and jB=0 when w=0. 

%% References:
% P. Kunder, "power system stability and control", 1994.

%% 
function Ybus = YbusCalc(ListLine,varargin) 

GridType = 'AC';        % Default is AC
for n = 1:length(varargin)
    if(strcmpi(varargin{n},'GridType'))
        GridType = varargin{n+1};   % 1-AC, 2-DC
    end
end

% Get the data
fb = ListLine(:,1);             % From bus number
tb = ListLine(:,2);             % To bus number

r = ListLine(:,3);              % Resistance,  R
x = ListLine(:,4);              % Inductance,  wL
b = ListLine(:,5);              % Capacitance, wC
g = ListLine(:,6);              % Conductance, G

T = ListLine(:,7);              % Turns ratio, T

% Calculate y
if strcmpi(GridType,'AC')
zs = r + 1j*x;                   
yp = g + 1j*b;       	% g and b can be "inf" without causing problems
z  = zs + 1./yp;      	% Total impedance of that branch
y  = 1./z;           	% Total admittance of that branch
elseif strcmpi(GridType,'DC')
    y = 1./r;           % For DC grid power flow, only the resistance is considerred, as w0 is 0 which means jX=0 and jB=0.
else
    error(['Error: GridType.'])
end

% Get the number
n_bus = max(max(fb),max(tb));    % Number of buses
n_branch = length(fb);           % Number of branches, including self branches

% Initialise YBus
Ybus = zeros(n_bus,n_bus);        

for k=1:n_branch     
    if(fb(k) ~= tb(k))    
        % Formation of the Off Diagonal Elements...
        Ybus(fb(k),tb(k)) = Ybus(fb(k),tb(k)) - y(k)/T(k);      % Mutual admittance is negative
        Ybus(tb(k),fb(k)) = Ybus(fb(k),tb(k));
        
        % Formation of the Diagonal Elements...
        Ybus(fb(k),fb(k)) = Ybus(fb(k),fb(k)) + y(k)/(T(k)^2);  % Self admittance is positive
        Ybus(tb(k),tb(k)) = Ybus(tb(k),tb(k)) + y(k);
    else
        % Formation of the Diagonal Elements...
        Ybus(fb(k),tb(k)) = Ybus(fb(k),tb(k)) + y(k);
    end 
end
 
end