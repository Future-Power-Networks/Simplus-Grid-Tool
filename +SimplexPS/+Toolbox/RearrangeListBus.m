% This function rearranges the netlist data of buses (i.e., power flow
% settings).

% Author(s): Yitong Li

function [UpdateBus,N_Bus] = RearrangeListBus(UserData)

% Load data
ListBus = xlsread(UserData,'Bus');

% Re-order the bus sequence and area
[i1] = unique(ListBus(:,11));
N_Area = length(i1);
if N_Area ~= max(i1)
    error(['Error: The total area number is different from the maximum area index.'])
end
for i2 = 1:N_Area
    UpdateBusArea{i2} = ListBus(find(ListBus(:,11)==i1(i2)),:);
    UpdateBusArea{i2} = sortrows(UpdateBusArea{i2},1);
    i3 = unique(UpdateBusArea{i2}(:,12));
    if length(i3)~=1
        error(['Error: In an area, the area type has to be the same for all buses.']);
    end
end

% Error check
[N_Bus,ColumnMax_Bus] = size(ListBus);
if N_Bus~=length(ListBus(:,3))
    error(['Error: The total bus number is different from the maximum bus index.'])
end
[n1,n2] = mode(ListBus(:,1));
if n2~=1
    error(['Error: Bus ' num2str(n1) ' is defined multiple times.'])
end
if (ColumnMax_Bus>12)
    error(['Error: Bus data overflow.']) 
end
for i4 = 1:N_Area
    ListBusType = UpdateBusArea{i4}(:,2);
    IndexBusSlack = find(ListBusType == 1);
    if (isempty(IndexBusSlack))
        error(['Error: The system has no slack bus in area ' num2str(i4) '.']);
    elseif IndexBusSlack ~= 1
        error(['Error: The first bus in each area has to be the slack bus.']);
    elseif (length(IndexBusSlack) > 1)
        error(['Error: The system has more than one slack bus in area ' num2str(i4) '.']); 
    end
end

% Output bus data
UpdateBus = sortrows(ListBus,1);

end