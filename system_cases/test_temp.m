Ym=GmDSS_Cell{2}(1:2,1:2);
Ym0=evalfr(Ym,1i*2*pi*0.1)
%Ym0_norm=norm(Ym0)
Ym0_svd=svd(Ym0)

Zm0=inv(Ym0)
Zm0_svd=svd(Zm0)