% This function checks if the model dimension obtained from the descritptor
% state space model matches the dimension that obtained from strings.

% Author(s): Yitong Li

%%
function [lx,lu,ly] = ObjCheckDim(Obj)
    % Calculate the dimension from dss model
    [~,ModelDss] = Obj.GetDSS(Obj);
    [lx,lu,ly] = SimplusGT.DssGetDim(ModelDss);
    
    % Calculate the dimension from strings
    [StateString,InputString,OutputString] = Obj.GetString(Obj);
    lx_ = length(StateString);
    lu_ = length(InputString);
    ly_ = length(OutputString);
    
    % Compare
    if (lx_~=lx)
        error(['Error: state dimiension mismathch']); 
    end
    if (lu_~=lu)
        error(['Error: input dimiension mismathch']);
    end
    if (ly_~=ly)
        error(['Error: output dimiension mismathch']); 
    end
end