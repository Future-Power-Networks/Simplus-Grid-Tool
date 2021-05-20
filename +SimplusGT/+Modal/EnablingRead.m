function [StatePFEnable, BodeEnable, Layer12Enable, Layer3Enable] = EnablingRead(filename)
%this function is to read enabling settings of PF analysis.
%Author: Yue Zhu
Enable = xlsread(filename,3);
if Enable == 0
    error('Please enable at least one function of participation analysis.')
end
StatePFEnable = Enable(1);
BodeEnable = Enable(2);
Layer12Enable = Enable(3);
Layer3Enable = Enable(4);

end