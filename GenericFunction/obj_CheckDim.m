% This function checks if the model dimension obtained from the descritptor
% state space model matches the dimension that obtained from strings.

% Author(s): Yitong Li

%%
function [lx,lu,ly] = obj_CheckDim(obj)
    % Calculate the dimension from dss model
    [~,ModelDSS] = obj.ReadDSS(obj);
    [lx,lu,ly] = dss_GetDim(ModelDSS);
    
    % Calculate the dimension from strings
    [StateString,InputString,OutputString] = obj.ReadString(obj);
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