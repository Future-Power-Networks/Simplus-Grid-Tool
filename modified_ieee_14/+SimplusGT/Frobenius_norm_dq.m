% Frobinus norm of a 2*2 d-q matrix
% Author: Yue Zhu
% Norm_F = sqrt(trace(A*A^H)): A^H: conjugate transpose of A. 

function norm = Frobenius_norm_dq(a)

norm = sqrt( a.dd*conj(a.dd) + a.dq*conj(a.dq) + a.qd*conj(a.qd) + a.qq*conj(a.qq) );

end