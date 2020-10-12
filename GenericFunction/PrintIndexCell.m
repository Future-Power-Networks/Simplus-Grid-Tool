% This function prints the elements in a cell array.

% Author(s): Yitong Li

function [Index_n] = PrintIndexCell(StringCell,cmax,Index_1)
    if isempty(StringCell)
        Index_n = Index_1;
    else
        fprintf('        ')
        CountPrintStr = 0;
        ls = length(StringCell);
        Index = Index_1-1;
        for k = 1:ls
            Index = Index+1;
            fprintf(['(' num2str(Index) ') ']);
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
        Index_n = Index;
    end
end