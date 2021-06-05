% This function calculates the "steady" nodal admittance matrix.

% Author(s): Yunjie Gu, Yitong Li

%% Notes:
%
% There are two different forms of branch, and the parallel form is used
% because it fits better for the power system network line and load
% structure.
% Series form:
%              ---             ---C---
%   FromBus---|a:1|---R---L---|       |---ToBus
%              ---             ---G---
% Parallel form
%              ---     -----R---L-----
%   FromBus---|a:1|---|               |---ToBus
%              ---    |    ---C---    |
%                      ---|       |---
%                          ---G---
% where (a:1) is the turn ratio of transformer.
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

% Get the data
FB = ListLine(:,1);             % From bus number
TB = ListLine(:,2);             % To bus number

R = ListLine(:,3);              % Resistance, R
X = ListLine(:,4);              % Inductance, wL
B = ListLine(:,5);              % Capacitance, wC
G = ListLine(:,6);              % Conductance, G

T = ListLine(:,7);              % Turns ratio, T

AreaTypeLine = ListLine(:,8);       % AC or DC type

% Get number
N_Bus = max(max(FB),max(TB));    % Number of buses
N_Branch = length(FB);

% Calculate y
BranchConnection = 2;      % 1-series; 2-parallel
for i = 1:N_Branch
    if AreaTypeLine(i) == 1
        Zs = R(i) + 1i*X(i);           
        Yp = G(i) + 1i*B(i);
        if BranchConnection == 1      % Series format
            Z  = Zs + 1/Yp;
            Y(i)  = 1/Z;
        elseif BranchConnection == 2  % Parallel format
            Y(i) = Yp + 1/Zs;
        else
            error(['Error: branch format.']);
        end
            
    elseif AreaTypeLine(i) == 2
        if R(i)==0 && isinf(G(i))
            error(['Error: Branch ' num2str(FB(i)) '-' num2str(TB(i)) ' is a DC branch, whose resistance can NOT be zero.'])
        end
        % For DC grid power flow, only the resistance is considerred, as w0 is 0 which means jX=0 and jB=0.
        if BranchConnection == 1
            Y(i) = 1/(1/G(i) + R(i)); 	
        elseif BranchConnection == 2
            Y(i) = G(i) + 1/R(i);
        end
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