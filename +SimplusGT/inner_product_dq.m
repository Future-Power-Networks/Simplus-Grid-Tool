% Frobinus inner product of a and b^H
% ^H:conjugate transpose of A. 
% Author: Yue Zhu


function product = inner_product_dq(a,b)

product = a.dd*b.dd + a.qd*b.dq + a.dq*b.qd + a.qq*b.qq;

end