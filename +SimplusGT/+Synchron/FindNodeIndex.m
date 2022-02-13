% Find the index of the first voltage bus/node
NumVbus1st = 1;

% Find the index of the first floating bus/node
NumFbus1st = N_Bus+1;  % Default: no floating bus
for i = 1:N_Bus
    if ApparatusSourceType(i) == 3
        NumFbus1st = i;
        break;
    end
end

% Find the index of the first current bus/node
NumIbus1st = N_Bus+1;   % Default: no current bus
for i = 1:N_Bus
    if ApparatusSourceType(i) == 2
        NumIbus1st = i;
        break;
    end
end

% Update bus index
if NumIbus1st>N_Bus
    if NumFbus1st <= N_Bus
        NumIbus1st = NumFbus1st;     % Update Ibus_1st
    end
end