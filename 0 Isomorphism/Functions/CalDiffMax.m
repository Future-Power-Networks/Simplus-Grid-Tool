function [DiffMax,Pos] = CalDiffMax(A)

[m,n] = size(A);

if m~=1 && n~=1
    error(['Error: Input is not an array.']);
end

imax = max(m,n);

DiffMax = 0;
Pos = [1,1];
for i1 = 1:imax
    a1 = A(i1);
    for i2 = i1:imax
        a2 = A(i2);
        Diff1 = abs(a2-a1);
        if Diff1>DiffMax
            DiffMax = Diff1;
            Pos = [i1,i2];
        end
    end
end

end