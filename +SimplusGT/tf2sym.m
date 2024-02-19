% This function converts a transfer function model to its symbolic form.

% Author(s): Yunjie Gu, Yitong Li

%%
% Notes:

% This function is based on tf which reported overflow when layout node>=5,
% and was therefore discarded

% This is much faster than ss2sym.m

% The input system can be a SISO transfer function or a MIMO transfer
% function matrix.

%%
function Gsym = tf2sym(Gtf)

    % Symbolic form of Laplace operator "s"
    s = sym('s');
    
    % Get the size of Gtf
    [rmax,cmax] = size(Gtf);
    
    % The input system is empty
    if (rmax <= 0) || (cmax <= 0)
        error(['Error: Empty input.']);
        
    % The input system is SISO
    elseif (rmax == 1) && (cmax == 1)
        [Num,Den] = tfdata(Gtf);
        Gsym = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s);
        
    % The input system is MIMO
    else
        for r = 1:rmax
            for c = 1:cmax
                Gtf_ = Gtf(r,c);
             	[Num,Den] = tfdata(Gtf_);
                Gsym(r,c) = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s);
            end
        end
    end
    
end