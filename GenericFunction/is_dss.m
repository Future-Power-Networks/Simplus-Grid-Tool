% This function checks if the system is in descriptor state space form.

% Author(s): Yitong Li

function Flag = is_dss(Gdss)

    A = Gdss.A;
    E = Gdss.E;
    
    if isempty(E) && (~isempty(A)) 
        Flag = 0;
    else
        Flag = 1;
    end
    
end