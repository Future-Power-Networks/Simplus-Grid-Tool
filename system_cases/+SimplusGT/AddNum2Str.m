% This function is used to add number index to a char string.

% Author(s): Yitong Li

function UpdateStr = AddNum2Str(Str,Num)

[r_Num,c_Num] = size(Num);
if r_Num~=1 && c_Num~=1
    error(['Error: The number has to be a scalar or a vector.']);
end
[r_Str,~] = size(Str);
if r_Str~=1
    error(['Error: The string has to be a 1*n cell.'])
end

for k = 1:length(Str)
    for i = 1:length(Num)
        Str{k} = [Str{k},num2str(Num(i))];
        
        % Notes:
        % This is for interlink apparatus, which is connected to multiple
        % buses. Some of its variables are not aligned to a specific bus,
        % for example, its dc-link voltage, between ac and dc sides.
        if i~=length(Num)
            Str{k} = [Str{k},'-'];
        end
    end
end

UpdateStr = Str;

end