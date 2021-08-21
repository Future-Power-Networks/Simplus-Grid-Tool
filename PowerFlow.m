% Program for Gauss - Seidel Load Flow Analysis
% Praviraj P G, MTech I Year, EE Dept., IIT Roorkee, India, Email :pravirajpg@gmail.com
% Modified by Y Gu. Adding self-conductance for passive load

% Assumption, Bus 1 is considered as Slack bus.
function [V,I,Av,P,Q,Vm]=PowerFlow(busdata,linedata)

ybus = YbusCalc(linedata);  % Calling program "YbusCalc.m" to get Y-Bus. 
bus = busdata(:,1);         % Bus number.
type= busdata(:,2);         % Type of Bus 1-Slack, 2-PV, 3-PQ
V   = busdata(:,3);         % Initial Bus Voltages.
th  = busdata(:,4);         % Initial Bus Voltage Angles.
GenMW   = busdata(:,5);     % PGi, Real Power injected into the buses.
GenMVAR = busdata(:,6);     % QGi, Reactive Power injected into the buses.
LoadMW  = busdata(:,7);     % PLi, Real Power Drawn from the buses.
LoadMVAR= busdata(:,8);     % QLi, Reactive Power Drawn from the buses.
Qmin = busdata(:,9);        % Minimum Reactive Power Limit
Qmax = busdata(:,10);       % Maximum Reactive Power Limit
nbus = max(bus);            % To get no. of buses
P = GenMW - LoadMW;         % Pi = PGi - PLi, Real Power at the buses.
Q = GenMVAR - LoadMVAR;     % Qi = QGi - QLi, Reactive Power at the buses.
Vprev = V;
toler = 1;                  % Tolerence.
iteration = 0;              % iteration starting
tolermax  = 1e-9;
iteramax  = 1e4;
%tolermax  = 1e-12;
%iteramax  = 1e6;
while (toler>tolermax) && (iteration<iteramax)
    for i = 2:nbus
        sumyv = 0;
        for k = 1:nbus
            if i ~= k
                sumyv = sumyv + ybus(i,k)* V(k);  % Vk * Yik
            end
        end
        if type(i) == 2             % Computing Qi for PV bus
            Q(i) = -imag(conj(V(i))*(sumyv + ybus(i,i)*V(i)));
            if (Q(i) > Qmax(i)) || (Q(i) < Qmin(i))  % Checking for Qi Violation.
                if Q(i) < Qmin(i)   % Whether violated the lower limit.
                    Q(i) = Qmin(i);
                else                % No, violated the upper limit.
                    Q(i) = Qmax(i);
                end
                type(i) = 3;        % If Violated, change PV bus to PQ bus.
            end
        end
        V(i) = (1/ybus(i,i))*((P(i)-1j*Q(i))/conj(V(i)) - sumyv); % Compute Bus Voltages.
        if type(i) == 2 % For PV Buses, Voltage Magnitude remains same, but Angle changes.
            V(i) = pol2rect(abs(Vprev(i)), angle(V(i)));
        end
    end
    iteration = iteration + 1;                  % Increment iteration count.
    I = ybus*V;
    S = I.*conj(V);
    N = length(V);
    tolerV = max(abs(abs(V) - abs(Vprev)));     % Calculate V tolerance.
    tolerP = max(abs(real(S(2:N))-P(2:N)));     % Calculate P tolerance, exclude the slack terminal
    tolerQ = max(abs(imag(S(2:N))+Q(2:N)));     % Calculate Q tolerance, exclude the slack terminal
    toler  = max([tolerV,tolerP,tolerQ]);       % Calculate total tolerance
    Vprev = V;  % Vprev is required for next iteration,  V(i) = pol2rect(abs(Vprev(i)), angle(V(i)));
end             % End of while loop / Iteration

Av = angle(V);  % Final Bus Voltage Angles in rad.
Vm = abs(V);    % Final Bus Voltage Amplitude.
S = V.*conj(I); % Final Appearant Power, Generator Convention.
P = real(S);    % Final Active Power, Generator Convention.
Q = imag(S);    % Final Reactive Power, Generator Convention.

end

function rect = pol2rect(r,o)   % r = magnitude, o = angle in radians.
rect = r*cos(o) + 1j*r*sin(o);  % rect = real + j*imag
end

