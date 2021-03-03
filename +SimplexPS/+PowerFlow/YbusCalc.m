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
function Ybus = YbusCalc(ListLine) 

% GridType = 'AC';        % Default is AC
% for n = 1:length(varargin)
%     if(strcmpi(varargin{n},'GridType'))
%         GridType = varargin{n+1};   % 1-AC, 2-DC
%     end
% end

% Get the data
FB = ListLine(:,1);             % From bus number
TB = ListLine(:,2);             % To bus number

R = ListLine(:,3);              % Resistance, R
X = ListLine(:,4);              % Inductance, wL
B = ListLine(:,5);              % Capacitance, wC
G = ListLine(:,6);              % Conductance, G

T = ListLine(:,7);              % Turns ratio, T

AreaType = ListLine(:,9);       % AC or DC type

% Get number
N_Bus = max(max(FB),max(TB));    % Number of buses
N_Branch = length(FB);

% Calculate y
for i = 1:N_Branch
    if AreaType(i) == 1
        Zs = R(i) + 1i*X(i);                   
        Yp = G(i) + 1i*B(i); 	% g and b can be "inf" without causing problems
        Z  = Zs + 1/Yp;      	% Total impedance of that branch
        Y(i)  = 1/Z;           	% Total admittance of that branch
    elseif AreaType(i) == 2
        if (FB(i) ~= TB(i)) && R(i)==0
            error(['Error: Branch ' num2str(FB(i)) '-' num2str(TB(i)) ' is a DC branch, whose resistance can NOT be zero.'])
        end
        Y(i) = 1/(1/G(i) + R(i)); 	% For DC grid power flow, only the resistance is considerred, as w0 is 0 which means jX=0 and jB=0.
    end
end

% Initialise YBus
Ybus = zeros(N_Bus,N_Bus);        

for k = 1:N_Branch     
    if(FB(k) ~= TB(k))    
        % Formation of the Off Diagonal Elements...
        Ybus(FB(k),TB(k)) = Ybus(FB(k),TB(k)) - Y(k)/T(k);      % Mutual admittance is negative
        Ybus(TB(k),FB(k)) = Ybus(FB(k),TB(k));
        
        % Formation of the Diagonal Elements...
        Ybus(FB(k),FB(k)) = Ybus(FB(k),FB(k)) + Y(k)/(T(k)^2);  % Self admittance is positive
        Ybus(TB(k),TB(k)) = Ybus(TB(k),TB(k)) + Y(k);
    else
        % Formation of the Diagonal Elements...
        Ybus(FB(k),TB(k)) = Ybus(FB(k),TB(k)) + Y(k);
    end 
end
 
end