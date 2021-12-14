% This function checks if a matrix is a squre matrix

% Author(s): Yitong Li

function Flag = is_squre(M)

[r,c] = size(M);

if r == c
    Flag = 1;
else
    Flag = 0;
end

end