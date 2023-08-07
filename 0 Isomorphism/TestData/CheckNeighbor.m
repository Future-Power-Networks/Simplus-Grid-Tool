% This function check the neighbour value of a element in a matrix

% Author(s): Yitong Li

function Output = CheckNeighbor(M,position,value,mode)

% Get the size of matrix M
[rmax,cmax] = size(M);

% Get the position of the target element
r = position(1);
c = position(2);

% Default
% Notes: If all neighbors equal to the value, the output = 1, else then the
% output = 0.
Output = 1;

if c~= cmax
    if M(r,c+1) ~= value
        Output = 0;
    end
end
if c~=1
    if M(r,c-1) ~= value
        Output = 0;
    end
end
if r~=rmax
    if M(r+1,c) ~= value
        Output = 0;
    end
end
if r~=1
    if M(r-1,c) ~= value
        Output = 0;
    end
end

if mode == 1
elseif mode == 2
    % Notes: 
    % For mode 1, we only check the 
    %     *
    % *  Ref  *
    %     *
    % But for mode 2, we check
    % *   *   *
    % *  Ref  *
    % *   *   *
    % i.e., checking 4 more positions.
    if r~=rmax && c~=cmax
        if M(r+1,c+1) ~= value
            Output = 0;
        end
    end
    if r~=rmax && c~=1
        if M(r+1,c-1) ~= value
            Output = 0;
        end
    end
    if r~=1 && c~=cmax
        if M(r-1,c+1) ~= value
            Output = 0;
        end
    end
    if r~=1 && c~=1
        if M(r-1,c-1) ~= value
            Output = 0;
        end
    end
else
    error(['Error: Error mode.'])
end

end