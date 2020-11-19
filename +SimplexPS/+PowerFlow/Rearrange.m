% This function Rearrange_PowerFlow

% Author(s): Yitong

function ListPowerFlow = Rearrange(PowerFlow)

lpf = length(PowerFlow);

for i = 1:lpf
    ListPowerFlow(i,:) = [i,PowerFlow{i}];
end

end