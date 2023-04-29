clear all
clc
close all

%%
% GFM analysis
Vd = sym('Vd','real');
Vq = sym('Vq','real');

Vp = Vd + 1i*Vq;
Vn = conj(Vp);

Z_GFM_pn = [1i*Vp,1i*Vp;
            -1i*Vn,-1i*Vn];

T = [1,1i;
     1,-1i];

 Z_GFM = T^(-1)*Z_GFM_pn*T
       
% GFL analysis
Id = sym('Id','real');
Iq = sym('Iq','real');

Ip = Id + 1i*Iq;
In = conj(Ip);

Y_GFL_pn = [Ip,-Ip;
            -In,In];
        
Y_GFL = T^(-1)*Y_GFL_pn*T

% Notes:
% Just be careful about the convention and the reference direction of
% currents.

%%
syms L w s

Z_L = [(s+1i*w)*L,  0;
       0,           (s-1i*w)*L];
Z_Ldq = inv(T)*Z_L*T;
Z_Ldq = simplify(Z_Ldq)

Y_Ldq = inv(Z_Ldq);
Y_Ldq = simplify(Y_Ldq)

%%
syms Z
Z_dq = [0,-Z;
        0,0];

Z_pn = T*Z_dq*inv(T);
Z_pn = simplify(Z_pn)
