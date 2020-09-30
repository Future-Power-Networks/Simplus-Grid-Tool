% This function will re-arrange the costormized data in the netlist

% Author(s): Yitong Li

%% Function
function [UpdateBus,UpdateLine,UpdateDevice,N_Bus,N_Branch,N_Device] = RearrangeNetlist(NetlistBus,NetlistLine,NetlistDevice) 

%% Re-arrange netlist bus
% Bus data for power flow analysis
[N_Bus,ColumnMax_Bus] = size(NetlistBus);
if (ColumnMax_Bus>10); error(['Error: bus data overflow.']); end
list_bus_type = NetlistBus(:,2);
index_slack = find(list_bus_type == 1);
if (isempty(index_slack)); error(['Error: no slack bus.']);
elseif (length(index_slack) > 1); error(['Error: more than one slack bus.']); 
end

% Re-arrange the bus data
NetlistBus = sortrows(NetlistBus,1);

% Output
UpdateBus = NetlistBus;

%% Re-arrange netlist line
% Organize data
[N_Branch,ColumnMax_Line] = size(NetlistLine); 

FB  = NetlistLine(:,1);   % From bus
TB  = NetlistLine(:,2);   % To bus
R   = NetlistLine(:,3);
X   = NetlistLine(:,4);
B   = NetlistLine(:,5);
G   = NetlistLine(:,6);

% Check data overflow
if (ColumnMax_Line>6); error(['Error: line data overflow.']); end

% Check number of bus
N_Bus_ = max(max(FB), max(TB) );
if (N_Bus_ ~= N_Bus); error(['Error: number of buses.']); end   

% Replace NaN by inf
netlist_line_NaN = isnan(NetlistLine);
[r,c] = find(netlist_line_NaN == 1);  	% Find the index of "inf"
NetlistLine(r,c) = inf;

% Check short-circuit and open-circuit
for i = 1:N_Branch
    if ( isinf(R(i)) || isinf(X(i)) || ((B(i)==0)&&(G(i)==0)) )
        error(['Error: branch' num2str(FB(i)) num2str(TB(i)) ' is open circuit']);
    elseif ( (R(i)==0) && (X(i)==0) && (isinf(B(i)) || isinf(G(i))) )
        error(['Error: branch' num2str(FB(i)) num2str(TB(i)) ' is short circuit']);
    elseif ((R(i)<0) || (X(i)<0) || (B(i)<0) || (G(i)<0) )
        error(['Error: negative line paramters']);
   	end
end

% Re-arrange the data
% Ensure "From Bus" <= "To Bus" is always valid
for i = 1:N_Branch
    if FB(i) > TB(i);
        [TB(i),FB(i)] = deal(FB(i),TB(i));
        [NetlistLine(i,2),NetlistLine(i,1)] = deal(NetlistLine(i,1),NetlistLine(i,2));
    end
end
% Re-order the branch sequence
NetlistLine = sortrows(NetlistLine,2);
NetlistLine = sortrows(NetlistLine,1);

% Output
UpdateLine = NetlistLine;

%% Re-arrange netlist device
% Data for devices connected to buses
[N_Device,ColumnMax_Device] = size(NetlistDevice);
if (ColumnMax_Device>9); error(['Error: device data overflow.']); end
if (N_Device ~= N_Bus); error(['Error: number of buses is different from number of devices.']); end

% Re-arrange
NetlistDevice = sortrows(NetlistDevice,1);

% Output
UpdateDevice = NetlistDevice;

end