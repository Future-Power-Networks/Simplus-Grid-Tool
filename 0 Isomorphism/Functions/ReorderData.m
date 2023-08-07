% This function re-orders the system apparatus data, so that it can follow
% this format: voltage node, current node, floating node

% Find pure empty bus, i.e., no apparatus and no passive load
Index_Ebus = [];
for k = 1:N_Bus
    if (ListBus(k,5)==0) && (ListBus(k,6)==0) && (ListBus(k,7)==0) && (ListBus(k,8)==0)
        Index_Ebus = [Index_Ebus,ListBus(k,1)];
    end
end

% Get the device source type:
for i = 1:N_Bus
    if DeviceType{i}==1
      	DeviceSourceType(i) = 1;    % Voltage node
    elseif DeviceType{i}==11
    	DeviceSourceType(i) = 2;    % Current node
    elseif DeviceType{i}==100
        DeviceSourceType(i) = 3;    % Floating node     
    else
     	error(['Error']);
    end
end

% Based on the device source type, we get the new order
Index_Vbus = find(DeviceSourceType == 1);
Index_Ibus = find(DeviceSourceType == 2);
Index_Fbus = find(DeviceSourceType == 3);
Order_Old2New = [Index_Vbus,Index_Ibus,Index_Fbus]; % Convert old order to new
Order_Old2New_NoFbus = [Index_Vbus,Index_Ibus];
for i = 1:N_Bus
    Order_New2Old(Order_Old2New(i)) = i;            % Convert new order back to old
end
for i = 1:(N_Bus-length(Index_Fbus))
    Order_New2Old_NoFbus_(Order_Old2New(i)) = i;
end
CounterFbus = 0;
for i = 1:length(Order_New2Old_NoFbus_)
    if Order_New2Old_NoFbus_(i) ~= 0
        Order_New2Old_NoFbus(i-CounterFbus) = Order_New2Old_NoFbus_(i);
    else
        CounterFbus = CounterFbus+1;
    end
end

% Existance of Node
if ~isempty(Index_Vbus); Exist_Vbus = 1; else; Exist_Vbus = 0; end
if ~isempty(Index_Ibus); Exist_Ibus = 1; else; Exist_Ibus = 0; end
if ~isempty(Index_Fbus); Exist_Fbus = 1; else; Exist_Fbus = 0; end

% Re-order device source tyoe
DeviceSourceType = DeviceSourceType(:,Order_Old2New);

% Re-order power flow
V = V(Order_Old2New,:);
I = I(Order_Old2New,:);

% Re-order nodal admittance matrix
Ybus = Ybus(Order_Old2New,Order_Old2New);

% Re-order device para
for i = 1:N_Bus
    DeviceTypeNew{i} = DeviceType{Order_Old2New(i)};
    DeviceParaNew{i} = DevicePara{Order_Old2New(i)};
end