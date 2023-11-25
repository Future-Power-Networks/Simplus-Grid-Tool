% Twist G to H by partially swapping the input and output

% Author(s): Yitong Li

%%
% Notes:
%
% The original system is
% [y1] = [G11 G12]*[u1]
% [y2]   [G21 G22] [u2]
% => The twisted system is
% [y1] = [H11 H12]*[u1]
% [u2]   [H21 H22] [y2]

%%
function [H] = MatrixTwist(G,index)

[G11,G12,G21,G22] = SimplusGT.MatrixPartition(G,index,index);

H11 = G11 - G12*inv(G22)*G21;
H12 = G12*inv(G22);
H21 = -inv(G22)*G21;
H22 = inv(G22);

H = [H11,H12;
     H21,H22];

end