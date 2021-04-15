% This function plots the Nyquist diagram of the symbolic system

% Author(s): Yitong Li

function nyquist_c(X,sbd,varargin)

[Half,~] = LoadVar(0,'Half',varargin);          % Default 0, i.e., plot bode rather than singular value

funcX = matlabFunction(X);

Xw = zeros(1,length(sbd));

for n = 1:length(sbd)
    try 
        Xw(n) = funcX(sbd(n));
    catch
        Xw(n) = funcX();
    end
end

W = length(Xw);

if (imag(sbd(1))*imag(sbd(W)))<0
wbdn = sbd(1:W/2);
wbdp = sbd(W/2+1:W);

Xwn  = Xw(1:W/2);
Xwp  = Xw(W/2+1:W);
Xwn = flip(Xwn);

Arg_wn = angle(Xwn);
Arg_wp = angle(Xwp);

for n = 1:W/2
    Xwn_x(n) = real(Xwn(n));
    Xwn_y(n) = imag(Xwn(n));
    Xwp_x(n) = real(Xwp(n));
    Xwp_y(n) = imag(Xwp(n));
end

plot(Xwn_x,Xwn_y); grid on; hold on;

if Half == 0
    plot(Xwp_x,Xwp_y); grid on; hold on;
end

else
    error(['Error: Please add both positive and negative frequency']);
end

end