% Author(s): Yitong Li

clear all
clc
close all

DataBus = importdata('Data300BusRawBus.txt');
DataBranch = importdata('Data300BusRawBranch.txt');

% DataBus
% 1       | 2   | 3   | 4   | 5        | 6 | 7     | 8  | 9  | 10 | 11 | 12 
% Bus No. | NaN | NaN | NaN | Bus Type | V | theta | PL | QL | Pg | Qg | NaN
% 13  | 14   | 15   |
% NaN | Qmax | Qmin |

[rmax,cmax] = size(DataBus);

for i = 1:rmax
    DataBus(i,7) = DataBus(i,7)/180*pi; % Convert degree to rad
    DataBus(i,8) = DataBus(i,8)/100;    % Convert to 100MVA
    DataBus(i,9) = DataBus(i,9)/100;    % Convert to 100MVA
    DataBus(i,10) = DataBus(i,10)/100;    % Convert to 100MVA
    DataBus(i,11) = DataBus(i,11)/100;    % Convert to 100MVA
    if DataBus(i,5) == 3
        DataBus(i,5) = 1;
    end
    if DataBus(i,5) == 0
        DataBus(i,5) = 3;
    end
end

NumApp = length(find(DataBus(:,5)==2));

% DataBranch
% 15
% TurnRatio

[rmax,cmax] = size(DataBranch);
for i = 1:rmax
    if DataBranch(i,15) == 0
        DataBranch(i,15) = 1;
    end
end

