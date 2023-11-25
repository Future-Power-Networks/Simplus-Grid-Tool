% Convert a nodal admittance matrix to hybrid
% admittance/impedance matrix H

% Author(s): Yitong Li

%%
% Notes:
%
% The original system is
% [I1] = [Y11 Y12]*[V1]
% [I2]   [Y21 Y22] [V2]
% => The hybrid system is
% [I1] = [H11 H12]*[V1]
% [V2]   [H21 H22] [I2]
%
% vbus_end is the index of the final voltage bus/node.

%%
function [H] = Y2H(Y,vbus_end)

[Y11,Y12,Y21,Y22] = SimplusGT.PartitionMatrix(Y,vbus_end,vbus_end);

H11 = Y11 - Y12*inv(Y22)*Y21;
H12 = Y12*inv(Y22);
H21 = -inv(Y22)*Y21;
H22 = inv(Y22);

H = [H11,H12;
     H21,H22];

end