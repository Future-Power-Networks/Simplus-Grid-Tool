function B = rgt_inv(A) %right inverse of a matrix

    B = SimplexPS.lft_inv(A.');
    B = B.';
    
end