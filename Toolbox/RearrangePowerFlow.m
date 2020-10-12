% This function Rearrange_PowerFlow

% Author(s): Yitong

function ListPowerFlow = RearrangePowerFlow(PowerFlow)

lpf = length(PowerFlow);

for i = 1:lpf
    ListPowerFlow(i,:) = [i,PowerFlow{i}];
end

end