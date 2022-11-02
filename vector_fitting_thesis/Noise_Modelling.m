Fs = 40000;
tx = 0:1/Fs:(n3-1)*1/Fs;
nx = 0:1:n3-1;
thetax=rand(1,8); %0~2pi????
wgnvd=sqrt(14.5)*randn(1,n3); %?????
wgnid=sqrt(0.02)*randn(1,n3);

nfx1=299.125*2*n3/fs+n3; %300Hz
nfx2=249.25*2*n3/fs+n3; %250Hz
nfx3=99.75*2*n3/fs+n3; %100Hz
nfx4=49.875*2*n3/fs+n3; %50Hz
nfx5=349.0*2*n3/fs+n3; %350Hz
nfx6=199.375*2*n3/fs+n3; %200Hz
nfx7=7852.625*2*n3/fs+n3;%7850
nfx8=8152.25*2*n3/fs+n3;%8150

vd_m300=sqrt(PSDvdabs(nfx1))*sin(2*pi*300*tx+thetax(1)*2*pi);   %300Hz
vd_m250=sqrt(PSDvdabs(nfx2))*sin(2*pi*250*tx+thetax(2)*2*pi);   %250Hz
vd_m100=sqrt(PSDvdabs(nfx3))*sin(2*pi*100*tx+thetax(3)*2*pi);   %100
vd_m50=sqrt(PSDvdabs(nfx4))*sin(2*pi*50*tx+thetax(4)*2*pi);    %50
vd_m350=sqrt(PSDvdabs(nfx5))*sin(2*pi*350*tx+thetax(5)*2*pi);   %350
vd_m200=sqrt(PSDvdabs(nfx6))*sin(2*pi*200*tx+thetax(6)*2*pi);   %200
vd_m7850=sqrt(PSDvdabs(nfx7))*sin(2*pi*7850*tx+thetax(7)*2*pi);  %7850
vd_m8150=sqrt(PSDvdabs(nfx8))*sin(2*pi*8150*tx+thetax(8)*2*pi);  %8150

id_m300=sqrt(PSDidabs(nfx1))*sin(2*pi*300*tx+thetax(1)*2*pi); %300Hz
id_m250=sqrt(PSDidabs(nfx2))*sin(2*pi*250*tx+thetax(3)*2*pi); %250Hz
id_m100=sqrt(PSDidabs(nfx3))*sin(2*pi*100*tx+thetax(3)*2*pi);   %100
id_m50=sqrt(PSDidabs(nfx4))*sin(2*pi*50*tx+thetax(4)*2*pi);    %50
id_m350=sqrt(PSDidabs(nfx5))*sin(2*pi*350*tx+thetax(5)*2*pi);   %350
id_m200=sqrt(PSDidabs(nfx6))*sin(2*pi*200*tx+thetax(6)*2*pi);   %200
id_m7850=sqrt(PSDidabs(nfx7))*sin(2*pi*7850*tx+thetax(7)*2*pi);  %7850
id_m8150=sqrt(PSDidabs(nfx8))*sin(2*pi*8150*tx+thetax(8)*2*pi);  %8150


vd_model = vd_m300+vd_m250+vd_m100+vd_m50+vd_m350+vd_m200+vd_m7850+vd_m8150+wgnvd;
id_model = id_m300+wgnid;%id_m250+id_m100+id_m50+id_m350+id_m200+id_m7850+id_m8150+wgnid;
N = length(tx);
ACFvd_m = xcorr(vd_model)/N;
ACFid_m = xcorr(id_model)/N;
CCFvdid_m=xcorr(vd_model,id_model)/N;
Pxid = fftshift(fft(ACFid_m));
Pxvd = fftshift(fft(ACFvd_m));

figure(10)
clf;
subplot(1,3,1);
plot((-N+1:N-1),ACFvd_m,'r');
% hold on;
% plot(-n3+1:n3-1,ACFvd,'b');
subplot(1,3,2);
plot((-N+1:N-1),ACFid_m,'r');
% hold on;
% plot(-n3+1:n3-1,ACFid,'b');
subplot(1,3,3);
plot((-N+1:N-1),CCFvdid_m,'r');
% hold on;
% plot(-n3+1:n3-1,CCFvdid,'b');
% 
% figure(12)
% plot((-N+1:N-1)*Fs/(2*N-1),abs(Px2)/Fs),'r';
% grid on;
% title('PSD');

figure(13)
plot(tx,vd_model,'r')
hold on;
%title('modeled noise')
plot(tx,vdx0,'b')
%title('Real noise')
%axis([1,1.05,-10,10])