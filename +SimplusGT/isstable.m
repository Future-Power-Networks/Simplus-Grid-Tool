% Author(s): Yunjie Gu
% Copyright, Simplus, Ltd.

function [y,u] = isstable(G,varargin)
    
    e = 1e-6;
    if ~isempty(varargin)
        e = varargin{1};
    end
    p = vpasolve(1/G);
    
    u = [];
    y = 1;
    for n =1:length(p)
        if(real(p(n))>e)
            u = [u;p(n)]; %#ok<AGROW>
            y = 0;
        end
    end
end