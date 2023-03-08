% This function prints the elements in a array.

% Author(s): Yitong Li

function PrintArray(Data,cmax)
    fprintf('        ')
    CountPrint = 0;
    ls = length(Data);
    
    for k = 1:ls
        fprintf(num2str(Data(k)))
        CountPrint = CountPrint + 1;
        if k ~= ls
            fprintf(', ');
            if CountPrint >= cmax
                fprintf('\n')
                CountPrint = 0;
                fprintf('        ')
            end
        else
            fprintf('\n');
        end
    end
end