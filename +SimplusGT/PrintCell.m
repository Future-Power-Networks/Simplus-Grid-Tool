% This function prints the elements in a cell array.

% Author(s): Yitong Li

function PrintCell(StringCell,cmax)
    fprintf('        ')
    CountPrintStr = 0;
    ls = length(StringCell);
    for k = 1:ls
        fprintf(StringCell{k})
        CountPrintStr = CountPrintStr + 1;
        if k ~= ls
            fprintf(', ');
            if CountPrintStr >= cmax
                fprintf('\n')
                CountPrintStr = 0;
                fprintf('        ')
            end
        else
            fprintf('\n');
        end
    end
end