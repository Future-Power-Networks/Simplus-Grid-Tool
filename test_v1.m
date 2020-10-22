
close all
% GH = tf([-1,0],[1])
% D = 1+GH;
% H = 1/D;
% pole_H = pole(H)
% 
% fn = fn+1;
% figure(fn)
% nyquist(GH)
% title('GH')
% 
% fn = fn+1;
% figure(fn)
% nyquist(D)
% title('D')
% 
% fn = fn+1;
% figure(fn)
% nyquist(H)
% title('H')

% s = sym('s');
% % s = tf('s');
% 

% s = sym('s');
% GH = -s;
% D = 1+GH;
% H = 1/D;

% omega_p = logspace(-5,5,1000)*2*pi;
% omega_pn = [-flip(omega_p),omega_p];
% 
% fn = fn+1;
% figure(fn)
% bodeplot(1)
% title('Mine')

t1 = ss([1],[1],[1],[1]);
t2 = t1;
t3 = t1;
t4 = t1;

Gcell{1,1} = t1;
Gcell{1,2} = t2;
Gcell{2,1} = t3;
Gcell{2,2} = t4;

t5 = ss_Arrange(Gcell)