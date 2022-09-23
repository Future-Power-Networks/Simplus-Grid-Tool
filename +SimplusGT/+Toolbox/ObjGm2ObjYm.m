% Author(s): Yitong Li

function ObjYm = ObjGm2ObjYm(ObjGm,NumBus)

% Get the string of apparatus
[GmStateStr,GmInStr,GmOutStr] = ObjGm.GetString(ObjGm);

% Get the index of electrical ports, i.e., voltage port and current port.
VecPortV = [];
VecPortI = [];

for i = 1:NumBus
    [~,AcPortV] = SimplusGT.CellFind(GmInStr,['v_d',num2str(i)]);
    [~,DcPortV] = SimplusGT.CellFind(GmInStr, ['v',num2str(i)]);
    
    [~,AcPortI] = SimplusGT.CellFind(GmOutStr,['i_d',num2str(i)]);
    [~,DcPortI] = SimplusGT.CellFind(GmOutStr,['i',num2str(i)]);
    if ~isempty(AcPortV)
        VecPortV = [VecPortV,AcPortV,AcPortV+1];
        VecPortI = [VecPortI,AcPortI,AcPortI+1];
        CellPortV{i} = [AcPortV,AcPortV+1];
        CellPortI{i} = [AcPortI,AcPortI+1];
    elseif ~isempty(DcPortV)
        VecPortV = [VecPortV,DcPortV];
        VecPortI = [VecPortI,DcPortI];
        CellPortV{i} = DcPortV;
        CellPortI{i} = DcPortI;
    else
        error(['Error']);
    end
end

ObjYm = SimplusGT.ObjTruncate(ObjGm,VecPortI,VecPortV);

end

