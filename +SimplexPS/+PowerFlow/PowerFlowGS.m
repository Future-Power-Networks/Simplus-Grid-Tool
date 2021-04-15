% Program for Gauss - Seidel Power Flow Analysis

function [PowerFlow,Ybus,V,I,Ang0,P,Q,Vm]=PowerFlowGS(ListBus,ListLine,w0)

Ybus = SimplexPS.PowerFlow.YbusCalc(ListLine);      % Get nodal admittance matrix
list_number = ListBus(:,1);     % Bus number
n_bus = max(list_number);       % Total number of buses
list_type = ListBus(:,2);      	% Bus type: 1-Slack, 2-PV, 3-PQ

index_slack = find(list_type == 1);      % Index of slack bus

if (isempty(index_slack))
    error(['Error: no slack bus']);
elseif (length(index_slack) > 1)
    error(['Error: more than one slack bus']);
elseif (index_slack ~= 1)
    error(['Error: bus 1 is not slack bus']);
end
number_slack = list_number(index_slack);

V0   = ListBus(:,3);         % Initial bus voltages.
th0  = ListBus(:,4);         % Initial bus voltage angles.

PGi = ListBus(:,5);     % Active power injected into the buses, G-generation
QGi = ListBus(:,6);     % Reactive power injected into the buses, G-generation
PLi = ListBus(:,7);     % Active power drawn from the buses, L-load
QLi = ListBus(:,8);     % Reactive power drawn from the buses, L-load

Qmin = ListBus(:,9);        % Minimum Reactive Power Limit
Qmax = ListBus(:,10);       % Maximum Reactive Power Limit

P = PGi - PLi;  	% Net actove power at buses.
Q = QGi - QLi;      % Net reactive power at buses.

%V = pol2rect(V0,th0);    	% Convert voltages from polar form to rectangular form
V = V0;
Vprev = V;

tolerance = 1;           	% Initialize tolerence
iteration = 0;              % Initialize interaction count

tolerance_max  = 1e-9;
iteration_max  = 1e4;

while ((tolerance>tolerance_max) && (iteration<iteration_max))
    
    for i = 1:n_bus
        if i ~= number_slack    % Check if slack bus
            
            sum_yv = 0;     
            for k = 1:n_bus
                if i ~= k
                    sum_yv = sum_yv + Ybus(i,k)* V(k);  % Vk * Yik
                end
            end
        
            if list_type(i) == 2             % Computing Qi for PV bus
                Q(i) = -imag(conj(V(i))*(sum_yv + Ybus(i,i)*V(i)));     % Equation (6.91) in Kunder's book
                if (Q(i) > Qmax(i)) || (Q(i) < Qmin(i))  % Checking for Qi violation
                    if Q(i) < Qmin(i)
                        Q(i) = Qmin(i); % Set Qi to lower limit
                    else          
                        Q(i) = Qmax(i); % Set Qi to upper limit
                    end
                    list_type(i) = 3;        % If violated, change bus type from PV to PQ
                end
            end
        
            V(i) = (1/Ybus(i,i))*((P(i)-1j*Q(i))/conj(V(i)) - sum_yv);  % Compute bus boltage
                                                                        % Equation (6.90) in Kunder's book
            
            if list_type(i) == 2 
                V(i) = SimplexPS.pol2rect(abs(Vprev(i)), angle(V(i))); % For PV bus, voltage magnitude remains same, but angle changes.
            end
        
        end
    end
    
    iteration = iteration + 1;                  % Increment iteration count.
    
    I = Ybus*V;
    S = I.*conj(V);
    N = length(V);
    
    tolerV = max(abs(abs(V) - abs(Vprev)));     % Calculate V tolerance.
    tolerP = max(abs(real(S(2:N)) - P(2:N)));     % Calculate P tolerance, exclude the slack terminal
    tolerQ = max(abs(imag(S(2:N)) + Q(2:N)));     % Calculate Q tolerance, exclude the slack terminal
    tolerance  = max([tolerV,tolerP,tolerQ]);       % Calculate total tolerance
    Vprev = V;  % Vprev is required for next iteration,  V(i) = pol2rect(abs(Vprev(i)), angle(V(i)));
    
end             % End of while loop / Iteration

Ang0 = angle(V);  % Final Bus Voltage Angles in rad.
Vm = abs(V);    % Final Bus Voltage Amplitude.
S = V.*conj(I); % Final Appearant Power, Generator Convention.
P = real(S);    % Final Active Power, Generator Convention.
Q = imag(S);    % Final Reactive Power, Generator Convention.

for i = 1:n_bus
    % The negative signs make P and Q in load convention
    PowerFlow{i} = [-P(i) -Q(i) Vm(i) Ang0(i) w0];
end

end



