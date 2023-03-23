clear all
clc
close all

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