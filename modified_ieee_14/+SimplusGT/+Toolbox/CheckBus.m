% This function checks a bus type, area type, area number, based on
% ListBus.

% Author(s): Yitong Li

function [BusType,AreaNum,AreaType] = CheckBus(BusNum,ListBus)

Index = find(ListBus(:,1) == BusNum,1);

if isempty(Index)
    % If could not find this bus
    BusType  = [];
    AreaNum  = [];
    AreaType = [];
else
    % If finding the bus
    BusType  = ListBus(Index,2);
    AreaNum  = ListBus(Index,11);
    AreaType = ListBus(Index,12);
end

end