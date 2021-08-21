% created 2019.12.11 by Yunjie Gu
%
%                       --CCC--
% branch :  --RRR--LLL--       ----
%                       --GGG--
%
% netlist:   |  From |  To   |   R   |  wL   |  wC   |   G   |
%            |  Bus  |  Bus  |       |       |       |       |
% linedata = [  1       2        0       X       0      inf;
%               1       1        0       0       0       G;
%               2       2        0       0       0       0];

function Ybus = YbusCalc(linedata) 

fb = linedata(:,1);             % From bus number...
tb = linedata(:,2);             % To bus number...
r = linedata(:,3);              % Resistance,  R...
x = linedata(:,4);              % Inductance,  wL...
b = linedata(:,5);              % Capacitance, wC...
g = linedata(:,6);              % Conductance, G...

zs= r + 1j*x;                   
yp= g + 1j*b;                   
z = zs + 1./yp;
y = 1./z;

nbus = max(max(fb),max(tb));    % no. of buses...
nbranch = length(fb);           % no. of branches...
Ybus = zeros(nbus,nbus);        % Initialise YBus...

% Formation of the Off Diagonal Elements...
for k=1:nbranch
    if(fb(k) ~= tb(k))    
        Ybus(fb(k),tb(k)) = Ybus(fb(k),tb(k)) - y(k);
        Ybus(tb(k),fb(k)) = Ybus(tb(k),fb(k)) - y(k);
        Ybus(fb(k),fb(k)) = Ybus(fb(k),fb(k)) + y(k);
        Ybus(tb(k),tb(k)) = Ybus(tb(k),tb(k)) + y(k);
    else
        Ybus(fb(k),tb(k)) = Ybus(fb(k),tb(k)) + y(k);
    end 
end
 
end