% Author(s): Yitong Li

%% 
function Ybus = YbusCalc_s_sym(ListLine,W0,FrameFlag) 

s = sym('s');

% Get the data
FB = ListLine(:,1);             % From bus number
TB = ListLine(:,2);             % To bus number

R = ListLine(:,3);              % Resistance, R
X = ListLine(:,4);              % Inductance, wL
B = ListLine(:,5);              % Capacitance, wC
G = ListLine(:,6);              % Conductance, G

T = ListLine(:,7);              % Turns ratio, T

AreaType = ListLine(:,9);       % AC or DC type

XL = ListLine(:,8);           % Inductive load

% Get number
N_Bus = max(max(FB),max(TB));    % Number of buses
N_Branch = length(FB);

% Calculate y
for i = 1:N_Branch
    if strcmpi(FrameFlag,'albe')   % ab frame
        Zs = R(i) + X(i)/W0*(s);                   
        Yp = G(i) + B(i)/W0*(s); 	% g and b can be "inf" without causing problems
    elseif strcmp(FrameFlag,'dq')   % dq frame
      	Zs = R(i) + X(i)/W0*(s+1i*W0);                   
        Yp = G(i) + B(i)/W0*(s+1i*W0); 	% g and b can be "inf" without causing problems
    else
        error(['Error']);
    end
    Z  = Zs + 1/Yp;                 % Total impedance of that branch
    Y(i)  = 1/Z;                    % Total admittance of that branch
    if ~isinf(XL(i))
        if FB(i)~=TB(i)
            error(['Error: The mutual branch can not have XL normally.'])
        else
            Y(i) = Y(i) + 1/(s*XL(i)/W0);
        end
    end
end

% Initialise YBus
Ybus = zeros(N_Bus,N_Bus);

% Convert Ybus from double to symbolic
Ybus = sym(Ybus);

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