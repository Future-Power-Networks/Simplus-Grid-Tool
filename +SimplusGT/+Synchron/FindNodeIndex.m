% Find the index of the first voltage bus/node
NumVbus1st = 1;

% Find the index of the first floating bus/node
NumFbus1st = NumBus+1;  % Default: no floating bus
for i = 1:NumBus
    if ApparatusSourceType(i) == 3
        NumFbus1st = i;
        break;
    end
end

% Find the index of the first current bus/node
NumIbus1st = NumBus+1;   % Default: no current bus
for i = 1:NumBus
    if ApparatusSourceType(i) == 2
        NumIbus1st = i;
        break;
    end
end

% Update bus index
if NumIbus1st>NumBus
    if NumFbus1st <= NumBus
        NumIbus1st = NumFbus1st;     % Update Ibus_1st
    end
end