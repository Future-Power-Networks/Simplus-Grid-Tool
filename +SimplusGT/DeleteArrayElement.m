% This function deletes certain elements in an array A according to the
% index n.

% Author(s): Yitong Li

function [B] = DeleteArrayElement(A,n)

    % Sort n
    n = sort(n);
    
    % Initialize B
    B = A([1:(n(1)-1)]);
    
    % Get B
    [r,c] = size(A);
    for i = 1:length(n)
     	if (i)<length(n)
         	C = A((n(i)+1):(n(i+1)-1));
        else
            C = A((n(i)+1):end);
        end
        if (r>1) && (c==1)
            B = [B;
                 C];
        elseif (r==1) && (c>1)
            B = [B,C];
        else
            error(['Error: The input A is not an array.']);
        end
    end

end