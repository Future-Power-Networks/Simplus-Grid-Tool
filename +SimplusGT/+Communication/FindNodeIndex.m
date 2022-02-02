% Find the index of the first voltage bus/node
n_Vbus_1st = 1;

% Find the index of the first floating bus/node
n_Fbus_1st = N_Bus+1;  % Default: no floating bus
for i = 1:N_Bus
    if ApparatusSourceType(i) == 3
        n_Fbus_1st = i;
        break;
    end
end

% Find the index of the first current bus/node
n_Ibus_1st = N_Bus+1;   % Default: no current bus
for i = 1:N_Bus
    if ApparatusSourceType(i) == 2
        n_Ibus_1st = i;
        break;
    end
end

% Update bus index
if n_Ibus_1st>N_Bus
    if n_Fbus_1st <= N_Bus
        n_Ibus_1st = n_Fbus_1st;     % Update Ibus_1st
    end
end