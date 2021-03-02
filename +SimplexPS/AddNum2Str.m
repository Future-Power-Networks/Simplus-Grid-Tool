% This function is used to add number index to a char string.

% Author(s): Yitong Li

function UpdateStr = AddNum2Str(Str,Num)

[r_Num,c_Num] = size(Num);
if r_Num~=1 && c_Num~=1
    error(['The number has to be a scalar or a vector.']);
end
[r_Str,~] = size(Str);
if r_Str~=1
    error(['The string has to be a 1*n cell.'])
end

for k = 1:length(Str)
    for i = 1:length(Num)
        Str{k} = [Str{k},num2str(Num(i))];
        if i~=length(Num)
            Str{k} = [Str{k},'-'];
        end
    end
end

UpdateStr = Str;

end