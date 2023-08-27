% Author(s): Yitong Li

clear all
clc
close all

DataBusAll = importdata('IEEE118BusRawData.txt');
DataBranch = importdata('IEEE118BranchRawData.txt');
DataBus = DataBusAll.data;

% DataBus
% 1       | 3        | 4 | 5     | 6  | 7  | 8  | 9  | 12   | 13
% Bus No. | Bus Type | V | theta | PL | QL | Pg | Qg | Qmax | Qmin

[rmax,cmax] = size(DataBus);

for i = 1:rmax
    DataBus(i,5) = DataBus(i,5) - DataBus(69,5);    % Convert to slack base
end

for i = 1:rmax
    DataBus(i,5) = DataBus(i,5)/180*pi; % Convert degree to rad
    DataBus(i,6) = DataBus(i,6)/100;    % Convert to 100MVA
    DataBus(i,7) = DataBus(i,7)/100;    % Convert to 100MVA
    DataBus(i,8) = DataBus(i,8)/100;    % Convert to 100MVA
    DataBus(i,9) = DataBus(i,9)/100;    % Convert to 100MVA
    DataBus(i,12) = DataBus(i,12)/100;    % Convert to 100MVA
    DataBus(i,13) = DataBus(i,13)/100;    % Convert to 100MVA
    if DataBus(i,3) == 3
        DataBus(i,3) = 1;
    end
    if DataBus(i,3) == 0
        DataBus(i,3) = 3;
    end
end

NumApp = length(find(DataBus(:,3)==2));

% DataBranch
% 15
% TurnRatio

[rmax,cmax] = size(DataBranch);
for i = 1:rmax
    if DataBranch(i,15) == 0
        DataBranch(i,15) = 1;
    end
end

