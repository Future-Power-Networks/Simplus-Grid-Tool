% This function converts a nodal admittance matrix to hybrid
% admittance/impedance matrix

% Author(s): Yitong Li

%%
% Notes:
%
% The original system is
% [I1] = [Y11 Y12]*[V1]
% [I2]   [Y21 Y22] [V2]
% => The hybrid system is
% [I1] = [G11 G12]*[V1]
% [V2]   [G21 G22] [I2]
%
% Ibus_1st is the index of the first current bus/node.

%%
function [G] = HybridMatrixYZ(Y,Ibus_1st)

[Y11,Y12,Y21,Y22] = SimplusGT.PartitionMatrix(Y,Ibus_1st-1,Ibus_1st-1);

G11 = Y11 - Y12*inv(Y22)*Y21;
G12 = Y12*inv(Y22);
G21 = -inv(Y22)*Y21;
G22 = inv(Y22);

G = [G11,G12;
     G21,G22];

end