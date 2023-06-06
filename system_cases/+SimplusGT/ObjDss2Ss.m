% Convert dss model object to ss model object

% Author(s): Yitong Li

function ObjSs = ObjDss2Ss(ObjDss)
    [~,Gdss] = ObjDss.GetDSS(ObjDss);
    [StateStrDss,InputStrDss,OutputStrDss] = ObjDss.GetString(ObjDss);
    
    E = Gdss.E;
    
    Tor = 1e-15;
    Counter = 1;
    for i = 1:length(E)
        if E(i,i) >= Tor
            StateStrSs{Counter} = StateStrDss{i};
            Counter = Counter+1;
        end
    end
    Gss = SimplusGT.dss2ss(Gdss);
    
    ObjSs =SimplusGT.Class.ModelBase;
    ObjSs.SetSS(ObjSs,Gss);
    ObjSs.SetString(ObjSs,StateStrSs,InputStrDss,OutputStrDss);
    
end