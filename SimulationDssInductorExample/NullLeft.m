function Z = NullLeft(A)

    [m,n] = size(A);
    % Orthonormal basis 
    [U,S,V] = svd(A);
    U
    if isempty(A)
        Z = U;
    else
        r = rank(A);
        Z = U(:,r+1:n);
    end

end