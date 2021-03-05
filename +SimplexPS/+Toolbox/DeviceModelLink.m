% This function links the device models

% Author(s): Yitong Li, Yunjie Gu

function Gobj = DeviceModelLink(GmObj)

    % Create a new obj
    Gobj = SimplexPS.Class.ModelBase;
    Gobj.SetDSS(Gobj,dss([],[],[],[],[]));
    Gobj.SetString(Gobj,{},{},{});
    
    % Append
    for n = 1:length(GmObj)
        Gobj = SimplexPS.ObjAppend(Gobj,GmObj{n});
    end
    
end