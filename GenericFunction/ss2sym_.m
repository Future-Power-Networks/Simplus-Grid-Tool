% convert a state-space model to a symbolic transfer function

% Author(s): Yunjie Gu

%%
% this is very time consuming and should be used with care
% 2020-02-09 accelarated by diagnolizing A, much faster but still slow
% the accelarated version is acceptable and used to replace tf2sym to solve
% overflow problem of tf

%%
function G = ss2sym_(S)

s = sym('s','real');
A = S.A;
B = S.B;
C = S.C;
D = S.D;

G = C*inv(s*eye(length(A)) - A)*B + D;
 
end