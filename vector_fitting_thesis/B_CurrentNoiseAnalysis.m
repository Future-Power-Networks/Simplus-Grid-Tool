%% data load
clear;
noise_samp = importdata('P1200W_1.csv');
fs=20e3; %20 kHz
T=30;
Idn=noise_samp.data(:,3);
Iqn=noise_samp.data(:,4);
t=(0:1/fs:(T-1/fs));
N=T*fs;

%% autocorrelation and PSD
Idn_ac=Idn-mean(Idn);
Rd=xcorr(Idn_ac)/N;
IdPSD=fftshift(fft(Rd)/fs);
IdPSD_log=10*log10(abs(IdPSD));%change the unit from W to dBm
IdPSD_abs=abs(IdPSD);

figure(1);
clf;
subplot(2,1,1)
plot(t,Idn);
%hold on;
%plot(t,Iqn);
%legend('id','iq');
ylabel('current/A')
xlabel('time/s')
title('waveform of Id: 30 s');
grid on;
subplot(2,1,2)
txx=0:1/fs:(0.2-1/fs);
plot(txx, Idn(1:length(txx)));
%hold on;
%plot(txx, Iqn(1:length(txx)));
title('waveform of Id (zoomed in)');
ylabel('current/A')
xlabel('time/s')
grid on;

figure(2);
clf;
subplot(2,1,1)
plot(-N+1:N-1,Rd);
title('Noise Autocorrelation')
grid on;


subplot(2,1,2)
plot((0:N-1)/(2*N-1)*fs,IdPSD_abs(N:end));
title('Id noise PSD - W')
xlabel('Frequency / Hz');
ylabel('Noise PSD W/Hz');
grid on;

fAxis=(0:N-1)/(2*N-1)*fs;
%% Noise Modelling
j=1;
clear fm; % fm records the frequency points that is going to be modelled.
for i=N:length(IdPSD_abs)
    if IdPSD_abs(i)>=0.002
        if j>1 && abs(fAxis(i-N+1)-fm(j-1).freq_Hz)<0.5 % if the frequency difference with the last point is less than 0.5 Hz, choose the point with larger power
            if IdPSD_abs(i)>fm(j-1).Power % if larger, cover the previous
                fm(j-1).freq_Hz=fAxis(i-N+1);
                fm(j-1).PointCount=i;
                fm(j-1).Power=IdPSD_abs(i);
            else % smaller or equal, abandon this data point
            end
        else % new data point
        	fm(j).freq_Hz=fAxis(i-N+1);
            fm(j).PointCount=i;
            fm(j).Power=IdPSD_abs(i);
            j=j+1;
        end
    end
end
Nfm=length(fm);
Theta_i=rand(1,Nfm)*2*pi;
Noise_M=0;
sigma2=Rd(N);
for i=1:Nfm
    fm(i).Id_amp = sqrt(fm(i).Power /7.5);
    fm(i).Id_M= fm(i).Id_amp * sin( 2*pi*fm(i).freq_Hz*t+Theta_i(i) ); % A*sin(2*pi*f*t+theta)
    Noise_M=Noise_M+fm(i).Id_M;
    sigma2=sigma2-0.5*(fm(i).Id_amp)^2;
    fm(i).Power_real=fm(i).Power /7.5;
end
%Noise_M=fm(1).Id_M/7.5;
Id_WGNm=sqrt(sigma2)*randn(1,N); % white Gaussion Noise modelling
Noise_M=Noise_M+Id_WGNm;

figure(3);
clf;
%subplot(2,1,1);
plot(t,Idn);
%axis([20,20.2,-3.8,-3]);
grid on;
hold on;
%subplot(2,1,2);
plot(t,Noise_M+mean(Idn));
title('Modelled Noise of Id');
axis([20,20.2,-3.8,-3])
grid on;
legend('measured noise', 'modelled noise')
ylabel('current/A')
xlabel('time/s')

%% PSD plot of modelled noise
Rd_M=xcorr(Noise_M)/N;
IdPSD_M=fftshift(fft(Rd_M)/fs);
IdPSD_M_log=10*log10(abs(IdPSD_M));%change the unit from W to dBm
IdPSD_M_abs=abs(IdPSD_M);
figure(4);
clf;
subplot(3,1,1)
plot(-N+1:N-1,Rd_M);
title('Model: Autocorrelation of Id')
sigma=Rd_M(N-1); % white noise.
subplot(3,1,2)
plot((0:N-1)/(2*N-1)*fs,IdPSD_M_abs(N:end));
title('Model: Id noise PSD - W')
xlabel('Frequency / Hz');
ylabel('PSD W/Hz');
grid on;
subplot(3,1,3)
fAxis=(0:N-1)/(2*N-1)*fs;
plot(fAxis,IdPSD_M_log(N:end));
title('Model: Id noise PSD - Bd')
xlabel('Frequency / Hz');
ylabel('PSD Bd/Hz');
grid on;

save Noise_M.mat Noise_M
save Fmodelling.mat fm
save sigma2.mat sigma2
save Idn.mat Idn
