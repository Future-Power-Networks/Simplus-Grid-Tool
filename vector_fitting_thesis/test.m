clear
t=0:0.001:1;
%y=sin(2*pi*10*t);
y(1)=1;
for i=2:length(t)
    if mod(i,500)==0
        y(i)=y(i-1)*(-1);
    else
        y(i)=y(i-1);
    end
end
plot(t,y);

