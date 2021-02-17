% This function convert the load power flow to the self branch impedance.

% Author(s): Yitong Li

function [UpdateListBus,UpdateListLine,UpdatePowerFlow] = Load2SelfBranch(ListBus,ListLine,DeviceType,PowerFlow)

%% Update "ListBus"
BusIndex  = ListBus(:,1);
BusType = ListBus(:,2);
PG      = ListBus(:,5);
QG      = ListBus(:,6);
PL      = ListBus(:,7);
QL      = ListBus(:,8);
N_Bus = max(BusIndex);

UpdateListBus = ListBus;

% Set PL and QL to zero
UpdateListBus(:,7) = zeros(size(ListBus(:,7)));
UpdateListBus(:,7) = zeros(size(ListBus(:,7)));

%% Update "PowerFlow"
UpdatePowerFlow = PowerFlow;
for i = 1:N_Bus
    % P and Q are in load convention
    P(i) = PowerFlow{i}(1);
    Q(i) = PowerFlow{i}(2);
    V(i) = PowerFlow{i}(3);
    
    % Update "PowerFlow"
    % Remove PL and QL from the power flow
    UpdatePowerFlow{i}(1) = P(i) - PL(i);
    UpdatePowerFlow{i}(2) = Q(i) - QL(i);
end

%% Error check
for i = 1:N_Bus
    if (PG(i)==0) && (QG(i)==0) && (DeviceType{i}~=100) && (BusType(i)~=1)
        error(['Error: Bus ' num2str(i) ' should be a slack bus (power flow) or floating bus (device) because PGi=0 and QGi=0.']);
    end
    if PL(i) < 0
        error(['Error: Passive load at bus ' num2str(i) ' can not generate active power, i.e., PLi can not be less than 0.']);
    end
end

%% Update ListLine
FB  = ListLine(:,1);   % From bus
TB  = ListLine(:,2);   % To bus
N_Branch = length(FB);

% Initialize "ListLine" for inductive load
UpdateListLine = [ListLine,inf([N_Branch,1],'double')]; % Set all XL to inf defaultly
for i = 1:N_Bus
    
    % Assume the load is parallel RL or RC
    GL(i) = PL(i)/V(i)^2;
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
    if isinf(GL(i)) || XL(i)==0 || isinf(BL(i))
        error(['Error: The passive load at bus ' num2str(i) 'is short-circuit.']);
    end
    
    if GL(i)~=0 || ~isinf(XL(i)) || BL(i)~=0
        % if the branch is not open-circuit
        n = SimplexPS.Toolbox.FindBranch(UpdateListLine,i,i);
        if ~isempty(n)
            % n is NOT empty, which means ListLine has this branch
            UpdateListLine(n,5) = UpdateListLine(n,5) + BL(i);
            UpdateListLine(n,6) = UpdateListLine(n,6) + GL(i);
            UpdateListLine(n,8) = XL(i);
        else
            % n is empty, which means ListLine does not have this branch.
            % In this case, a new branch should be added
            UpdateListLine = [UpdateListLine;
                              i,i,0,0,BL(i),GL(i),1,XL(i)];
        end
    else
        % Do not need to care about the open-circuit
    end
end

% Re-order the branch sequence
UpdateListLine = sortrows(UpdateListLine,2);
UpdateListLine = sortrows(UpdateListLine,1);

end