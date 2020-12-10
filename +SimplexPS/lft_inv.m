function B = lft_inv(A) %left inverse of a matrix

    [n,m] = size(A);

    % QR = A, where Q is orthogonal, R is upper triangular (QR decomposition)  
    % R = [ P ] ---   m*m
    %     [ O ] --- (n-m)*m
    % S = [ P^(-1) O ]
    % SR = I
    % SQ^(-1) * A = SQ^(-1) * QR = I

    [Q,R] = qr(A);

    if rank(R) < m
        error('the matrix is not invertable');
    end

    P = R(1:m,:);
    S = [P^(-1) , zeros(m,n-m)];
    B = S * Q^(-1);

end