% This function rearranges the netlist data of buses (i.e., power flow
% settings).

% Author(s): Yitong Li

function [UpdateBus,N_Bus] = RearrangeListBusStruct(inputData)

% Load data
ListBus=[];
for i = 1:size(inputData.buses,1)
    bus=inputData.buses(i);
    ListBus(i,1) = bus.busNo;
    ListBus(i,2) = bus.busType;
    ListBus(i,3) = bus.VSp;
    ListBus(i,4) = bus.theta;
    ListBus(i,5) = bus.PGi;
    ListBus(i,6) = bus.QGi;
    ListBus(i,7) = bus.PLi;
    ListBus(i,8) = bus.QLi;
    ListBus(i,9) = bus.Qmin;
    ListBus(i,10) = bus.Qmax;
    ListBus(i,11) = bus.areaNo;
    ListBus(i,12) = bus.ACDC;
end

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
    if i3 == 2
        i4 = find(ListBus(:,2)==2,1);
        if ~isempty(i4)
            error(['Error: The bus type in dc network can not be 2!']);
        end
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
for i5 = 1:N_Bus
    if ListBus(i5,12) == 2
        if ListBus(i5,2) == 2
            error(['Error: Bus ' num2str(ListBus(i5,1)) ' is a dc bus, whose type can not be 2.'])
        end
    end
end

% Output bus data
UpdateBus = sortrows(ListBus,1);

end