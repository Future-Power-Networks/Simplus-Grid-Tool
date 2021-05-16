% convert a transfer function model to symbolic
% this is much faster than ss2sym.m

% Author(s): Yunjie Gu

%%
% 2020-02-09 this function is based on tf which reported overflow when 
% layout node>=5, and was therefore discarded

function G = tf2sym(S)

    s = sym('s');
    
    [Num,Den] = tfdata(S);
    G = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s);
    
end