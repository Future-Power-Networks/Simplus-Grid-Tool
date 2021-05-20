% Program for Newton-Raphson Load Flow Analysis..

function [PowerFlow] = PowerFlowNR(ListBus,ListLine,w0)

Y = SimplusGT.PowerFlow.YbusCalc(ListLine); 	% Calling ybusppg.m to get Y-Bus Matrix..
busd = ListBus;      	% Calling busdatas..

bus = busd(:,1);      	% Bus Number..
nbus = length(bus);   	% Number of buses

type = busd(:,2);      	% Type of Bus 1-Slack, 2-PV, 3-PQ..
V = busd(:,3);        	% Specified Voltage..
del = busd(:,4);     	% Voltage Angle..
Pg = busd(:,5);     	% PGi..
Qg = busd(:,6);     	% QGi..
Pl = busd(:,7);         % PLi..
Ql = busd(:,8);         % QLi..
Qmin = busd(:,9);       % Minimum Reactive Power Limit..
Qmax = busd(:,10);      % Maximum Reactive Power Limit..

P = Pg - Pl;         	% Pi = PGi - PLi..
Q = Qg - Ql;          	% Qi = QGi - QLi..
Psp = P;              	% P Specified..
Qsp = Q;             	% Q Specified..
G = real(Y);          	% Conductance matrix..
B = imag(Y);         	% Susceptance matrix..

pv = find(type == 2 | type == 1);   % PV Buses..
pq = find(type == 3);               % PQ Buses..
npv = length(pv);                   % No. of PV buses..
npq = length(pq);                   % No. of PQ buses..

Tol = 1;  
Iter = 1;
while (Tol > 1e-5)   % Iteration starting..
    
    P = zeros(nbus,1);
    Q = zeros(nbus,1);
    % Calculate P and Q
    for i = 1:nbus
        for k = 1:nbus
            P(i) = P(i) + V(i)* V(k)*(G(i,k)*cos(del(i)-del(k)) + B(i,k)*sin(del(i)-del(k)));
            Q(i) = Q(i) + V(i)* V(k)*(G(i,k)*sin(del(i)-del(k)) - B(i,k)*cos(del(i)-del(k)));
        end
    end

    % Checking Q-limit violations..
    if Iter <= 7 && Iter > 2    % Only checked up to 7th iterations..
        for n = 2:nbus
            if type(n) == 2
                QG = Q(n)+Ql(n);
                if QG < Qmin(n)
                    V(n) = V(n) + 0.01;
                elseif QG > Qmax(n)
                    V(n) = V(n) - 0.01;
                end
            end
         end
    end
    
    % Calculate change from specified value
    dPa = Psp-P;
    dQa = Qsp-Q;
    k = 1;
    dQ = zeros(npq,1);
    for i = 1:nbus
        if type(i) == 3
            dQ(k,1) = dQa(i);
            k = k+1;
        end
    end
    dP = dPa(2:nbus);
    M = [dP; dQ];       % Mismatch Vector
    
    % Jacobian
    % J1 - Derivative of Real Power Injections with Angles..
    J1 = zeros(nbus-1,nbus-1);
    for i = 1:(nbus-1)
        m = i+1;
        for k = 1:(nbus-1)
            n = k+1;
            if n == m
                for n = 1:nbus
                    J1(i,k) = J1(i,k) + V(m)* V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
                end
                J1(i,k) = J1(i,k) - V(m)^2*B(m,m);
            else
                J1(i,k) = V(m)* V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
            end
        end
    end
    
    % J2 - Derivative of Real Power Injections with V..
    J2 = zeros(nbus-1,npq);
    for i = 1:(nbus-1)
        m = i+1;
        for k = 1:npq
            n = pq(k);
            if n == m
                for n = 1:nbus
                    J2(i,k) = J2(i,k) + V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                J2(i,k) = J2(i,k) + V(m)*G(m,m);
            else
                J2(i,k) = V(m)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
    
    % J3 - Derivative of Reactive Power Injections with Angles..
    J3 = zeros(npq,nbus-1);
    for i = 1:npq
        m = pq(i);
        for k = 1:(nbus-1)
            n = k+1;
            if n == m
                for n = 1:nbus
                    J3(i,k) = J3(i,k) + V(m)* V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                J3(i,k) = J3(i,k) - V(m)^2*G(m,m);
            else
                J3(i,k) = V(m)* V(n)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
    
    % J4 - Derivative of Reactive Power Injections with V..
    J4 = zeros(npq,npq);
    for i = 1:npq
        m = pq(i);
        for k = 1:npq
            n = pq(k);
            if n == m
                for n = 1:nbus
                    J4(i,k) = J4(i,k) + V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
                end
                J4(i,k) = J4(i,k) - V(m)*B(m,m);
            else
                J4(i,k) = V(m)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
            end
        end
    end
    
    J = [J1 J2; J3 J4];     % Jacobian Matrix..
    
    X = inv(J)*M;           % Correction Vector
    dTh = X(1:nbus-1);      % Change in Voltage Angle..
    dV = X(nbus:end);       % Change in Voltage Magnitude..
    
    % Updating State Vectors..
    del(2:nbus) = dTh + del(2:nbus);    % Voltage Angle..
    k = 1;
    for i = 2:nbus
        if type(i) == 3
            V(i) = dV(k) + V(i);        % Voltage Magnitude..
            k = k+1;
        end
    end
    
    Iter = Iter + 1;
    Tol = max(abs(M));                  % Tolerance..
    
end

% [Pi,Qi] = LoadFlow(nbus,V,del,Sbase,Y,ListLine,ListBus);              % Calling Loadflow.m..
Sbase = 1;               	% Base, which is set to 1 currently
[Pi,Qi] = SimplusGT.PowerFlow.LoadFlow(nbus,V,del,Sbase,Y,ListLine,ListBus);              % Calling Loadflow.m..

for i = 1:nbus
    PowerFlow{i} = [-Pi(i) -Qi(i) V(i) del(i) w0];
end

end