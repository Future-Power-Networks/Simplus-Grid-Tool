% This function will re-arrange the costormized data in the netlist

% Author(s): Yitong Li

%% Function
function [UpdateBus,UpdateLine,UpdateDevice,N_Bus,N_Branch,N_Device] = RearrangeNetlist(ListBus,ListLine,ListDevice) 

%% Re-arrange netlist bus
% Bus data for power flow analysis
[N_Bus,ColumnMax_Bus] = size(ListBus);

% Re-order the bus sequence
ListBus = sortrows(ListBus,1);

% Output
UpdateBus = ListBus;

% Error check
if (ColumnMax_Bus>10)
    error(['Error: Bus data overflow.']) 
end
list_bus_type = ListBus(:,2);
index_slack = find(list_bus_type == 1);
if (isempty(index_slack))
    error(['Error: The system has no slack bus.']);
elseif index_slack ~= 1
    error(['Error: Bus 1 has to be the slack bus']);
elseif (length(index_slack) > 1)
    error(['Error: The system can not have more than one slack bus.']); 
end

%% Re-arrange netlist line
% Organize data
[N_Branch,ColumnMax_Line] = size(ListLine); 

FB  = ListLine(:,1);   % From bus
TB  = ListLine(:,2);   % To bus
Rbr = ListLine(:,3);
Xbr = ListLine(:,4);
Bbr = ListLine(:,5);
Gbr = ListLine(:,6);
Tbr = ListLine(:,7);

% Check data overflow
if (ColumnMax_Line>7)
    error(['Error: Line data overflow.']); 
end

% Check number of bus
N_Bus_FT = max(max(FB), max(TB) );
if (N_Bus_FT ~= N_Bus)
    error(['Error: Numbers of buses derived from "PowerFlow" and "NetworkLine" are different.'])
end   

% Replace NaN by inf
netlist_line_NaN = isnan(ListLine);
[r,c] = find(netlist_line_NaN == 1);  	% Find the index of "inf"
ListLine(r,c) = inf;

% Check short-circuit and open-circuit
for i = 1:N_Branch
    if ( isinf(Rbr(i)) || isinf(Xbr(i)) || ((Bbr(i)==0)&&(Gbr(i)==0)) )
    	error(['Error: Branch' num2str(FB(i)) num2str(TB(i)) ' is open circuit']);
    end
    if ( (Rbr(i)==0) && (Xbr(i)==0) && (isinf(Bbr(i)) || isinf(Gbr(i))) )
        error(['Error: Branch' num2str(FB(i)) num2str(TB(i)) ' is short circuit']);
    end
    if ((Rbr(i)<0) || (Xbr(i)<0) || (Bbr(i)<0) || (Gbr(i)<0) )
        error(['Error: Negative line paramters']);
    end
    if Tbr(i) <= 0
        error(['Error: Turns ratio can not be less than or equal to 0.']);
    end
end

% Re-arrange the data
% Ensure "From Bus" <= "To Bus" is always valid
for i = 1:N_Branch
    if FB(i) > TB(i)
        % Switch "From" and "To"
        [TB(i),FB(i)] = deal(FB(i),TB(i));
        [ListLine(i,2),ListLine(i,1)] = deal(ListLine(i,1),ListLine(i,2));
        
        % Switch the positions of line impedance and transformer
        Rbr(i) = Rbr(i)*Tbr(i)^2;
        ListLine(i,3) = ListLine(i,3)*ListLine(i,7)^2;
        Xbr(i) = Xbr(i)*Tbr(i)^2;
        ListLine(i,4) = ListLine(i,4)*ListLine(i,7)^2;
        
        % Change the turns ratio
        Tbr(i) = 1/Tbr(i);
        ListLine(i,7) = 1/ListLine(i,7);
    end
end

% Re-order the branch sequence
ListLine = sortrows(ListLine,2);
ListLine = sortrows(ListLine,1);

% Output
UpdateLine = ListLine;

%% Re-arrange netlist device
% Re-order the device sequence
ListDevice = sortrows(ListDevice,1);

% Output
UpdateDevice = ListDevice;

% Error check
[N_Device,ColumnMax_Device] = size(ListDevice);
if (ColumnMax_Device>50)
    error(['Error: Device data overflow.']); 
end
if (N_Device ~= N_Bus)
    error(['Error: Number of buses is different from number of devices.']); 
end

end