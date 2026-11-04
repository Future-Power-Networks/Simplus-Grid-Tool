function Flag = IsAcApparatusType(ApparatusType)
%ISACAPPARATUSTYPE True for active ac apparatus models.
%
% Types 0-89 are native white-box ac models. Types 4000 and 4001 are
% identified ac GFL and GFM admittance models, respectively. Either model
% can import dq- or sequence-domain impedance data.

Flag = ((ApparatusType >= 0) & (ApparatusType <= 89)) | ...
    ismember(ApparatusType,[4000,4001]);

end
