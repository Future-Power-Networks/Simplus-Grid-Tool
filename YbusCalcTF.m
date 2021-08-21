% created 2019.11.12 by Yunjie Gu
% based on YbusCalc.m
% transfer function Ybus
% input: branch from-to format
% branches with the same from-to are in paralell
% branches with from=to are self-admittances

function [Ybtf,Ybss] = YbusCalcTF(linedata,w) 

s = tf('s');

fb = linedata(:,1);             % From bus number...
tb = linedata(:,2);             % To bus number...
r = linedata(:,3);              % Resistance,  R...
x = linedata(:,4);              % Inductance,  wL...
b = linedata(:,5);              % Capacitance, wC...
g = linedata(:,6);              % Conductance, G...

nbus = max(max(fb),max(tb));        % no. of buses...
nbranch = length(fb);               % no. of branches...
Ybus = zeros(2*nbus,2*nbus)*s;      % Initialise YBus...
y = zeros(2,2,nbranch)*s;           % Initialise y branch...    

for n = 1:nbranch

    dp = s*b(n)/w + g(n);
    qp = b(n);
    ds = s*x(n)/w + r(n);
    qs = x(n);
    ps = ds^2 + qs^2;
    pp = dp^2 + qp^2;
    dt = ds*pp + dp;
    qt = qs*pp - qp;
    pt = 1 + pp*ps + 2*dp*ds - 2*qp*qs;
    
    if (isinf(r(n)) || isinf(x(n)) || ((g(n)==0) && (b(n)==0)))
        y(:,:,n) = zeros(2,2);
    elseif (isinf(g(n)) || isinf(b(n)))
        if ((r(n)==0) && (x(n)==0))
            error(['branch ' num2str(n) ' short circuit']);
        else            
            ys = ps^(-1)*[ds qs;-qs ds];
            ys = minreal(ys);
            y(:,:,n) = ys;
        end
    else
        yt = pt^(-1)*[dt qt;-qt dt];
        yt = minreal(yt);
        y(:,:,n) = yt;
    end
    
end

for k=1:nbranch
    if fb(k) ~= tb(k) %off diagonal
        Ybus((2*fb(k)-1):(2*fb(k)),(2*tb(k)-1):(2*tb(k))) = Ybus((2*fb(k)-1):(2*fb(k)),(2*tb(k)-1):(2*tb(k)))-y(:,:,k);
        Ybus((2*tb(k)-1):(2*tb(k)),(2*fb(k)-1):(2*fb(k))) = Ybus((2*tb(k)-1):(2*tb(k)),(2*fb(k)-1):(2*fb(k)))-y(:,:,k);
        Ybus((2*fb(k)-1):(2*fb(k)),(2*fb(k)-1):(2*fb(k))) = Ybus((2*fb(k)-1):(2*fb(k)),(2*fb(k)-1):(2*fb(k)))+y(:,:,k);
        Ybus((2*tb(k)-1):(2*tb(k)),(2*tb(k)-1):(2*tb(k))) = Ybus((2*tb(k)-1):(2*tb(k)),(2*tb(k)-1):(2*tb(k)))+y(:,:,k);
    else              %diagonal
        Ybus((2*fb(k)-1):(2*fb(k)),(2*tb(k)-1):(2*tb(k))) = Ybus((2*fb(k)-1):(2*fb(k)),(2*tb(k)-1):(2*tb(k)))+y(:,:,k);
    end 
end

Ybtf = Ybus;
Ybss = ss(Ybtf);

end
