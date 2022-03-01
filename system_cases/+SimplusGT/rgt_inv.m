function B = rgt_inv(A) %right inverse of a matrix

    B = SimplusGT.lft_inv(A.');
    B = B.';
    
end