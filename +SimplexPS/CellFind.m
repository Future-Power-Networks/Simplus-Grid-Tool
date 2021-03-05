% This functio achieves the "find" for cell type.

% Author(s): Yitong Li

function [r_Out,c_Out] = CellFind(Cell,a)

[r_Cell,c_Cell] = size(Cell);

r_Out = [];
c_Out = [];

for r = 1:r_Cell
    for c = 1:c_Cell
        a_ = Cell{r,c};
        if ischar(a_) && ischar(a)
            if strcmp(a_,a)
                r_Out = [r_Out,r];
                c_Out = [c_Out,c];
            end
        elseif isnumeric(a_) && isnumeric(a)
            if ~isempty(find(a_==a, 1))
              	r_Out = [r_Out,r];
                c_Out = [c_Out,c];
            end
        end
    end
end


end