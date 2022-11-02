
Fs=1000; %????

n=0:1/Fs:1;
xn=cos(2*pi*40*n)+3*cos(2*pi*100*n)+randn(size(n));

nfft=1024;

cxn=xcorr(xn,'unbiased'); %??????????

CXk=fft(cxn,nfft);

Pxx=abs(CXk);

index=0:round(nfft/2-1);

k=index*Fs/nfft;

plot_Pxx=10*log10(Pxx(index+1));

plot(k,Pxx(index+1));