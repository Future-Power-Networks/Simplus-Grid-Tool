% This function re-orders the system apparatus data, so that it can follow
% this format: voltage node, current node, floating node

function [ApparatusSourceType,ApparatusParaNew,IndexEbus,IndexVbus,IndexIbus,IndexFbus,OrderOld2New,Ybus,Ybus_,ExistVbus,ExistIbus,ExistFbus,V,I] = ReorderData(ListBus,ApparatusType,V,I,Ybus,Ybus_,Para)

[N_Bus,~] = size(ListBus);

% Find pure empty bus, i.e., no apparatus and no passive load
IndexEbus = [];
for k = 1:N_Bus
    if (ListBus(k,5)==0) && (ListBus(k,6)==0) && (ListBus(k,7)==0) && (ListBus(k,8)==0)
        IndexEbus = [IndexEbus,ListBus(k,1)];
    end
end

% Get the apparatus source type:
% Notes: Remember to include other kinds of apparatuses later
for i = 1:N_Bus
    if ApparatusType{i}== 0 || ApparatusType{i}== 1
      	ApparatusSourceType(i) = 1;    % Voltage node
    elseif ApparatusType{i}==11
    	ApparatusSourceType(i) = 2;    % Current node
    elseif ApparatusType{i}==100
        ApparatusSourceType(i) = 3;    % Floating node     
    else
     	error(['Error: The apparatus type ']);
    end
end

% Based on the device source type, we get the new order
IndexVbus = find(ApparatusSourceType == 1);
IndexIbus = find(ApparatusSourceType == 2);
IndexFbus = find(ApparatusSourceType == 3);
OrderOld2New = [IndexVbus,IndexIbus,IndexFbus]; % Convert old order to new
OrderOld2NewNoFbus = [IndexVbus,IndexIbus];
for i = 1:N_Bus
    OrderNew2Old(OrderOld2New(i)) = i;            % Convert new order back to old
end
for i = 1:(N_Bus-length(IndexFbus))
    OrderNew2OldNoFbus_(OrderOld2New(i)) = i;
end
CounterFbus = 0;
for i = 1:length(OrderNew2OldNoFbus_)
    if OrderNew2OldNoFbus_(i) ~= 0
        OrderNew2OldNoFbus(i-CounterFbus) = OrderNew2OldNoFbus_(i);
    else
        CounterFbus = CounterFbus+1;
    end
end

% Existance of Node
if ~isempty(IndexVbus); ExistVbus = 1; else; ExistVbus = 0; end
if ~isempty(IndexIbus); ExistIbus = 1; else; ExistIbus = 0; end
if ~isempty(IndexFbus); ExistFbus = 1; else; ExistFbus = 0; end

% Re-order device source tyoe
ApparatusSourceType = ApparatusSourceType(:,OrderOld2New);

% Re-order power flow
V = V(OrderOld2New,:);
I = I(OrderOld2New,:);

% Re-order nodal admittance matrix
Ybus = Ybus(OrderOld2New,OrderOld2New);
Ybus_ = Ybus_(OrderOld2New,OrderOld2New);

% Re-order device para
for i = 1:N_Bus
    ApparatusTypeNew{i} = ApparatusType{OrderOld2New(i)};
    ApparatusParaNew{i} = Para{OrderOld2New(i)};
end

end