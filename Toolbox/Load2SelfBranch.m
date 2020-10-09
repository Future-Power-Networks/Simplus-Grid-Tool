% This function convert the load power flow to the self branch impedance.

% Author(s): Yitong Li

function [UpdateListBus,UpdateListLine,UpdatePowerFlow] = Load2SelfBranch(ListBus,ListLine,DeviceType,PowerFlow)

%% Initialize Output
UpdateListBus = ListBus;
UpdateListLine = ListLine;
UpdatePowerFlow = PowerFlow;

%% Organize data
IndBus  = ListBus(:,1);
TypeBus = ListBus(:,2);
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
    if (PG(i)==0) && (QG(i)==0) && (DeviceType{i}~=100) && (TypeBus{i}~=1)
        error(['Error: Bus ' num2str(i) ' should be a slack bus (power flow) or floating bus (device) because PGi=0 and QGi=0.']);
    end
    if PL(i) < 0
        error(['Error: Passive load at bus ' num2str(i) ' can not generate active power, i.e., PLi can not be less than 0.']);
    end
end

%% Calculate the load value
for i = 1:N_Bus
    % Assume the load is parallel RL or RC
    RL(i) = V(i)^2/PL(i);
    if QL(i) > 0
        % RL load
        XL(i) = V(i)^2/QL(i);
        BL(i) = 0;
    elseif QL(i) <= 0
        % RC load
        XL(i) = inf;
        BL(i) = -QL(i)/V(i)^2;
    end
    
    % Error check
    if (RL(i) == 0) || (XL(i) == 0) || isinf(BL(i))
        error(['Error: The passive load at bus ' num2str(i) 'is short-circuit.']);
    end
end

% Update "ListLine"
UpdateListLine = [ListLine,inf([N_Branch,1],'double')]; % Set XL to inf defaultly
for i = 1:N_Branch
    if FB(i) == TB(i)
        UpdateListLine(i,5) = UpdateListLine(i,5)+BL(FB(i));    % Combine capacitive part to self branch
        UpdateListLine(i,6) = UpdateListLine(i,6)+1/RL(FB(i));  % Combine resistive part to self branch
        UpdateListLine(i,8) = XL(FB(i));                        % Update the inductive part
    end
end

end