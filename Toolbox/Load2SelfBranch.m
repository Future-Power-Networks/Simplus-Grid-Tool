% This function convert the load power flow to the self branch impedance.

% Author(s): Yitong Li

function [UpdateListBus,UpdateListLine,UpdatePowerFlow] = Load2SelfBranch(ListBus,ListLine,DeviceType,PowerFlow,LoadConnection)

%% Initialize Output
UpdateListBus = ListBus;
UpdateListLine = ListLine;
UpdatePowerFlow = PowerFlow;

%% Organize data
IndBus  = ListBus(:,1);
PG      = ListBus(:,5);
QG      = ListBus(:,6);
PL      = ListBus(:,7);
QL      = ListBus(:,8);
N_Bus = max(IndBus);

% Update "ListBus"
UpdateListBus(:,7) = zeros(size(ListBus(:,7)));
UpdateListBus(:,7) = zeros(size(ListBus(:,7)));

FB  = ListLine(:,1);   % From bus
TB  = ListLine(:,2);   % To bus
Rbr = ListLine(:,3);
Xbr = ListLine(:,4);
Bbr = ListLine(:,5);
Gbr = ListLine(:,6);
N_Branch = length(FB);

for i = 1:N_Bus
    % P and Q are in load convention
    P(i) = PowerFlow{i}(1);
    Q(i) = PowerFlow{i}(2);
    V(i) = PowerFlow{i}(3);
    
    % Update "PowerFlow"
    UpdatePowerFlow{i}(1) = P(i) - PL(i);
    UpdatePowerFlow{i}(2) = Q(i) - QL(i);
end

%% Error check
for i = 1:N_Bus
    if PG(i)==0 && QG(i)==0
        if DeviceType{i} ~= 100
            error(['Error: a floating bus is required according to the power flow analysis.']);
        end
    elseif PL(i) < 0
        error(['Error: passive load should not generate active power.']);
    elseif QL(i) < 0
        error(['Error: passive load should not generate reactive power.']);
    end
end

%% Calculate the load value
for i = 1:N_Bus
    
    % Series RL connection
    if LoadConnection == 1
        SL = PL(i)+j*QL(i);
        I = conj(SL/V(i));
        Z = V(i)/I;
        RL(i) = real(Z);
        XL(i) = imag(Z);
    % Parallel RL connection: This one is used
    elseif LoadConnection == 2
        RL(i) = V(i)^2/PL(i);
        XL(i) = V(i)^2/QL(i);
    else
        error(['Error']);
    end
    
    % Error check
    if (RL(i) == 0) || (XL(i) == 0)
        error(['Error']);
    end
end

% Update "ListLine"
UpdateListLine = [ListLine,zeros(N_Branch,1)];
for i = 1:N_Branch
    if FB(i) == TB(i)
        UpdateListLine(i,6) = UpdateListLine(i,6)+1/RL(FB(i));
        UpdateListLine(i,7) = XL(FB(i));
    end
end

end