function Z = NullRight(A)

    [m,n] = size(A);
    
    % Orthonormal basis 
    [U,S,V] = svd(A);
    V
    S
    if isempty(A)
        Z = V;
    else
        r = rank(A);
        Z = V(:,r+1:n);
    end

end