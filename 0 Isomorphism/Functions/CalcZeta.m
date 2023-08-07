% This function calculates zeta based on T and xi.

% Author(s): Yitong Li

function [zeta_m,w_min] = CalcZeta(T,xi)

s = sym('s');
w = sym('w','real');
T = subs(T,'s',1i*w);

w_min = [];
zeta_m = [];

for i = 1:length(xi)        
    
    zeta0{i} = xi(i)/s + 1/T;
        % Notes: zeta0{i} cannot be written as (xi(i) + s/T)/s, even though
        % it is equivalent to xi(i)/s + 1/T. Do not know why.
    zeta0{i} = subs(zeta0{i},'s',1i*w);

    zeta{i} = zeta0{i}*conj(zeta0{i});
    zeta{i} = simplify(zeta{i});

    zeta0_{i} = matlabFunction(zeta0{i});
    zeta_{i} = matlabFunction(zeta{i});

    w_min(i) = fminbnd(zeta_{i},-1e4,1e4);
    zeta_m(i) = abs(zeta0_{i}(w_min(i)));
    
end

end