% This function finds the index of a branch in ListLine

function [n] = FindBranch(ListLine,FB,TB)

FB_List = ListLine(:,1);
TB_List = ListLine(:,2);

N_Branch = length(FB_List);

n = [];
Count = 0;
for i = 1:N_Branch
    if FB == FB_List(i) && TB == TB_List(i)
        n = i;
        Count = Count+1;
    end
end

if Count >= 2
    error(['Error: The branch' num2str(FB) '-' num2str(TB) ' appears multiple times.'])
end

end