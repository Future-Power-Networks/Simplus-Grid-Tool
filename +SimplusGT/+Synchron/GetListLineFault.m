% This function converts a normal network line matrix "ListLine" to a fault
% network line matrix "ListLineFault" based on bus number of short-circuit
% fault.
%
% Author(s): Yitong Li

function [ListLineFault] = GetListLineFault(NfaultBus,ListLine,ListBus)

% Find fault bus
IndexFaultBus = find(ListBus(:,1) == NfaultBus);
if isempty(IndexFaultBus)
    error(['Error: The fault bus number does not exist in the system.'])
else
    FaultAreaType = ListBus(IndexFaultBus,12);
end

% Find fault branch
IndexBranchFrom = find(ListLine(:,1) == NfaultBus);
IndexBranchTo = find(ListLine(:,2) == NfaultBus);


[Index,~] = intersect(IndexBranchFrom,IndexBranchTo);

% Create a fault branch
% The fault can not be a pure short-circuit because the modeling will
% report an error.
%             | From | To | R | X | C | G | T | AreaType
ListLineAdd = [NfaultBus, NfaultBus, inf, inf, 0, 1e9, 1, FaultAreaType];

% Add the fault branch into the ListLine Matrix
if isempty(Index)
    ListLineFault = [ListLine;
                     ListLineAdd];
    ListLineFault = sortrows(ListLineFault,2);
    ListLineFault = sortrows(ListLineFault,1);
else
    ListLineFault = ListLine;
    ListLineFault(Index,:) = ListLineAdd;
end

end