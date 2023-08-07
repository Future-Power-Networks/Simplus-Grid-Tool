% This function prepares the new Ybus matrix for adding the inductors of
% SGs into the networks by D-Y conversion.

% Author(s): Yitong Li

function Ybus = PrepareConvertDY(Ybus,Ibus_1st,N_Bus,Ysg)

% Initialize Ybus
Ybus = blkdiag(Ybus,zeros(Ibus_1st-1));

for i = 1:Ibus_1st-1
    
    % Add bus i2 and add a branch between bus i1 and i2
    i1 = i;
    i2 = N_Bus+i;
    
    % Self branch
    Ybus(i1,i1) = Ybus(i1,i1) + Ysg{i1};
    Ybus(i2,i2) = Ybus(i2,i2) + Ysg{i1};
    
    % Mutual branch
    Ybus(i1,i2) = Ybus(i1,i2) - Ysg{i1};
    Ybus(i2,i1) = Ybus(i2,i1) - Ysg{i1};
    
end

% Switch bus position
Ybus = [Ybus(N_Bus+1:end,:);        % New voltage buses
        Ybus(Ibus_1st:N_Bus,:);     % Old current bus
        Ybus(1:Ibus_1st-1,:)];      % Old voltage bus (assume as zero current bus)
      % New voltage buses,  old current bus,       old voltage bus
Ybus = [Ybus(:,N_Bus+1:end),Ybus(:,Ibus_1st:N_Bus),Ybus(:,1:Ibus_1st-1)];

end