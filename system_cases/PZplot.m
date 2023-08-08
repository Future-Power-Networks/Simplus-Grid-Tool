Vs=1;
Zth=0.1+1i*1;
s= 10:-0.001:0.001;
for i=1:length(s)
    Zk(i) = s(i)-0.1*1i*s(i);
end
Vk = abs(Vs*(Zk./(Zth+Zk)));
Pk = abs((Vk).^2./Zk);
figure(2)
plot(Pk,(Vk))