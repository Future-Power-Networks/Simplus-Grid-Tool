% This function seperates the matrix into four sub matrices.

% Author(s): Yitong Li

function [M11,M12,M21,M22] = PartitionMatrix(M,r,c)

M11 = M(1:r,1:c);
M12 = M(1:r,(c+1):end);
M21 = M((r+1):end,1:c);
M22 = M((r+1):end,(c+1):end);

end