% Author(s): Yunjie Gu

function Y = sym2tf(X)

    [M,N]=size(X);
    Y = zeros(M,N)*tf(1,1);
    
    for m = 1:M
        for n = 1:N
            x = X(m,n);
            [xn,xd] = numden(x);
            xnp = sym2poly(xn);
            xdp = sym2poly(xd);
            %xnp = real(xnp);
            %xdp = real(xdp);
            y = tf(xnp,xdp);
            Y(m,n) = y;
        end
    end
    
    Y = minreal(Y);
end