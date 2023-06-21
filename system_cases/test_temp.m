Ym=GmDSS_Cell{2}(1:2,1:2);
Ym0=evalfr(Ym,-1i*2*pi*0.1);
%Ym0_norm=norm(Ym0)
%Ym0_svd=svd(Ym0)

Zm0=inv(Ym0);
Zm0_svd=svd(Zm0);

T=[1,1i; 1, -1i]/sqrt(2);
T*Zm0*T'

% 
% Ym_tf=tf(GmDSS_Cell{2}(1:2,1:2));
% 
% Zm_tf_pn = inv(T)*Zm_tf*T; % seq admittance, with FCE
% 
% evalfr(Zm_tf_pn, -1i*2*pi*0.1)


% Ym_tf_pp = Ym_tf_pn(1,1);
% Ym_tf_p = Ym_tf_pp - evalfr(Ym_tf_pp,1i*2*pi*60);
% evalfr(Ym_tf_pn, 1i*2*pi*60)

%Ym_tf_pn0 = evalfr(Ym_tf_pn(1,1), 1i*2*pi*0.1-1i*2*pi*60)


