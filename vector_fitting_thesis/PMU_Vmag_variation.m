%% data config
%PMU_str = importdata('UOS_PMU_01_12_21.csv');
tstart_hour = 10;
tend_hour = 18;
t_duration_sec=30;
n_start = tstart_hour*3600*50+1;
n_end = tend_hour*3600*50;
%n_end = (tstart_hour+(t_duration_sec/3600))*3600*50;
% fs=50;
% N = length(n_start:n_end);
% PMU_Va = PMU_str.data(n_start:n_end,3);
% PMU_Vb = PMU_str.data(n_start:n_end,5);
% PMU_Vc = PMU_str.data(n_start:n_end,7);
% t_s = (n_start:n_end)/fs;
% Vmag = sqrt(PMU_Va.^2+PMU_Vb.^2+PMU_Vc.^2);
% Freq = PMU_str.data(n_start:n_end,1);
%Vmag = Vmag + 5*cos(2*pi*100*t_s);

%% Probability Distrubution
clear V_diff;
for kx=n_start:(n_end-t_duration_sec*50)
    k_end= kx + t_duration_sec*50;
    Va=PMU_str.data(kx:k_end,3);
    Vb=PMU_str.data(kx:k_end,5);
    Vc=PMU_str.data(kx:k_end,7);
    Vmag=sqrt(Va.^2+Vb.^2+Vc.^2);
    V_diff_k = std(Vmag);
    V_diff(kx-n_start+1)=V_diff_k;
end

figure(1);
clf;
%subplot(1,2,1);
Bar_num=120;
[h,k] = hist(V_diff,Bar_num);
bar(k,h/trapz(k,h));

% %% Probability Distribution
% figure(1);
% clf;
% subplot(1,2,1);
% Bar_num=75;
% [h,k] = hist(Vmag,Bar_num);
% bar(k,h/trapz(k,h));
% title('Vmag probability distribution') ;
% xlabel('V');
% ylabel('probability');
% subplot(1,2,2);
% Bar_num=75;
% [h,k] = hist(Freq,Bar_num);
% bar(k,h/trapz(k,h));
% title('Frequency probability distribution') ;
% xlabel('Hz');
% ylabel('probability');
% 
% 
% %% Vmag across the time
% figure(2);
% clf;
% subplot(1,2,1);
% t_h=(n_start:n_end)/3600/50;
% plot(t_h,Vmag)
% title('Vmag') ;
% xlabel('hour');
% ylabel('V');
% subplot(1,2,2);
% plot(t_h,Freq)
% title('Frequency') ;
% xlabel('hour');
% ylabel('Hz');
% 
% %% Power density
% %figure(3);
% %periodogram(Vmag,rectwin(N),N,fs);
% 
% % Nfft=1024; % N-point fft
% % Vmag_AC=xcorr(Vmag,'unbiased');
% % Vmag_fft = fft(Vmag_AC,Nfft);
% % Vmag_fft_abs = abs(Vmag_fft);
% % index=round(Nfft/2-1);
% % f_point = index * fs / Nfft;
% % figure(3)
% % plot((0:n-1)*fs/(2*n-1),ACVmag_abs(n:end));
% % 
% figure(3);
% clf;
% subplot(1,2,1);
% Vmag_ = Vmag - mean(Vmag);
% Vmag_AC = xcorr(Vmag_,'biased');
% PSDVmag=fftshift(fft(Vmag_AC));
% ACVmag_log=10*log10(abs(PSDVmag));%change the unit from W to dBm
% ACVmag_abs=abs(PSDVmag);
% plot((0:N-1)/(2*N-1)*fs,ACVmag_log(N:end));
% ylim([-150 50])
% title('autocorrelation Vmag noise PSD')
% xlabel('Frequency / Hz');
% ylabel('PSD dB/Hz');
% grid on;
% 
% subplot(1,2,2);
% Freq_ = Freq - mean(Freq);
% Vmag_AC = xcorr(Freq_,'biased');
% PSDVmag=fftshift(fft(Vmag_AC));
% ACVmag_log=10*log10(abs(PSDVmag));%change the unit from W to dBm
% ACVmag_abs=abs(PSDVmag);
% plot((0:N-1)/(2*N-1)*fs,ACVmag_log(N:end));
% ylim([-150 50])
% title('autocorrelation grid frequency noise PSD')
% xlabel('Frequency / Hz');
% ylabel('PSD dB/Hz');
% grid on;
% 
% %Rxx = xcorr(x,'biased');
% %Rxxdft = abs(fftshift(fft(Rxx)));
% %freq = -Fs/2:Fs/length(Rxx):Fs/2-(Fs/length(Rxx));
% %plot(freq,Rxxdft);
% % xdft = fft(Vmag);
% % xdft = xdft(1:N/2+1);
% % psdx = (1/(Fs*N)) * abs(xdft).^2;
% % psdx(2:end-1) = 2*psdx(2:end-1);
% % freq = 0:Fs/N:Fs/2;
% % 
% % figure(3);
% % plot(freq,10*log10(psdx))
% % grid on
% % title('Periodogram Using FFT')
% % xlabel('Frequency (Hz)')
% % ylabel('Power/Frequency (dB/Hz)')
% 
