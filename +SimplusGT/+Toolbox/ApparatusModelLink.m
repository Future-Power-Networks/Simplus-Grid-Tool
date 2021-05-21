% This function links the apparatus models

% Author(s): Yitong Li, Yunjie Gu

function Gobj = ApparatusModelLink(GmObj)

    % Create a new obj
    Gobj = SimplusGT.Class.ModelBase;
    Gobj.SetDSS(Gobj,dss([],[],[],[],[]));
    Gobj.SetString(Gobj,{},{},{});
    
    % Append
    for n = 1:length(GmObj)
        Gobj = SimplusGT.ObjAppend(Gobj,GmObj{n});
    end
    
end