function MatrixNew = NormMatrixElement(MatrixOld,varargin)

% Set default value
DiagFlag = 1;

% Update
for n = 1:length(varargin)
    if(strcmpi(varargin{n},'DiagFlag'))
        DiagFlag = varargin{n+1};
    end
end

% Calculate the element norm
[rmax,cmax] = size(MatrixOld);
for r = 1:rmax
    for c = 1:cmax
        MatrixNew(r,c) = abs(MatrixOld(r,c));
    end
end

% Deal with the diagonal elements
if DiagFlag == 0
    % Set the diagonal elements of matrix to zero
    nmax = min(rmax,cmax);
    for n = 1:nmax
        MatrixNew(n,n) = 0;
    end
elseif DiagFlag == 1
    % Keep the diagonal elements of matrix
end

end