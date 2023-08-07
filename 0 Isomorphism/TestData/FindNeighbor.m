% This function finds certain value in its neighbor

% Author(s): Yitong Li

function [position_out] = FindNeighbor(M,position,value,mode)

% Get the size of the matrix
[rmax,cmax] = size(M);

% Get the position of the reference element
r = position(1);
c = position(2);

% Check mode first
if mode~=1 && mode~=2
    error(['Error: Error mode.']);
end

% Counter for output
counter = 0;

% -1 row
if r~=1 && c~=1 && mode==2
    if M(r-1,c-1) == value
        counter = counter + 1;
        position_out(counter,:) = [r-1,c-1]; 
    end
end
if r~=1
    if M(r-1,c) == value
        counter = counter + 1;
        position_out(counter,:) = [r-1,c];
    end
end
if r~=1 && c~=cmax && mode==2
    if M(r-1,c+1) == value
       counter = counter + 1;
       position_out(counter,:) = [r-1,c+1];
    end
end

% 0 row
if c~=1
    if M(r,c-1) == value
        counter = counter + 1;
        position_out(counter,:) = [r,c-1]; 
    end
end
if c~=cmax
    if M(r,c+1) == value
       counter = counter + 1;
       position_out(counter,:) = [r,c+1];
    end
end

% +1 row
if r~=rmax && c~=1 && mode==2
    if M(r+1,c-1) == value
        counter = counter + 1;
        position_out(counter,:) = [r+1,c-1]; 
    end
end
if r~=rmax
    if M(r+1,c) == value
       counter = counter + 1;
       position_out(counter,:) = [r+1,c];
    end
end
if r~=rmax && c~=cmax && mode==2
    if M(r+1,c+1) == value
       counter = counter + 1;
       position_out(counter,:) = [r+1,c+1];
    end
end

end