% This function find the mode in a cell.

% Author(s): Yitong

function [a,N_a] = CellMode(A)

[r_A,c_A] = size(A);

if r_A~=1 && c_A~=1
    error(['Error: The input A has to be a 1*n cell.'])
end

% Get the length
Length_A = length(A);

% Initialize the array
Array_A = [];

for i = 1:Length_A
    [r_Ai,c_Ai] = size(A{i});
    if r_Ai~=1 && c_Ai~=1
        error(['Error: The cell element of A has to be a 1*n array.']);
    end
    Length_Ai = length(A{i});
    for j = 1:Length_Ai
        Array_A = [Array_A,A{i}(j)];
    end
end

% Find the mode
[a,N_a] = mode(Array_A);

end