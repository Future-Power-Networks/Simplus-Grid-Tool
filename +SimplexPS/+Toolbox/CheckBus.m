% This function checks a bus type, area type, area number, based on
% ListBus.

function [BusType,AreaNum,AreaType] = CheckBus(BusNum,ListBus)

Index = find(ListBus(:,1) == BusNum,1);
BusType  = ListBus(Index,2);
AreaNum  = ListBus(Index,11);
AreaType = ListBus(Index,12);

end