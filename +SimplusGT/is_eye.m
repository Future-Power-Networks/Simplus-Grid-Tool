% This function checks if a matrix is an identity matrix.

% Author(s): Yitong Li

function Flag = is_eye(M)

Flag = true;

if (SimplusGT.is_squre(M) && isdiag(M))
    for i = 1:length(M)
        if M(i,i)~= 1
            Flag = 0;
        end
    end
else
    Flag = false;
end

end