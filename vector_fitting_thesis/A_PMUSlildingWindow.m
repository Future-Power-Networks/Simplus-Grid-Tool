%% data config
%clear;
%PMU_str = importdata('UOS_PMU_01_12_21.csv');
tstart_h = 10;
tend_h = 18;
t_length_sec=30;
n_start = tstart_h*3600*50+1;
n_end = tend_h*3600*50;
%n_end = (tstart_hour+(t_duration_sec/3600))*3600*50;
% fs=50;
% N = length(n_start:n_end);
PMU_Va = PMU_str.data(n_start:n_end,3);
PMU_Vb = PMU_str.data(n_start:n_end,5);
PMU_Vc = PMU_str.data(n_start:n_end,7);
% t_s = (n_start:n_end)/fs;
Vmag = sqrt(PMU_Va.^2+PMU_Vb.^2+PMU_Vc.^2)/sqrt(3);
Vmag_norm=Vmag./mean(Vmag);
% Freq = PMU_str.data(n_start:n_end,1);
%Vmag = Vmag + 5*cos(2*pi*100*t_s);
Vstd_mov=movstd(Vmag_norm,t_length_sec*50);
tx=tstart_h: 1/50/3600 : tend_h-1/50/3600;
std_fit=fitdist(Vstd_mov,'GeneralizedExtremeValue');
std_cdf=cdf(std_fit,0.001); % 0.1%

figure(1)
clf;
subplot(2,1,1)
plot(tx,Vmag_norm)
title('Normalised voltage magnitude from 10 am to 6 pm')
xlabel('time / hour')
ylabel('normalised voltage')
grid on

subplot(2,1,2)
plot(Vmag_norm(1:30),'linewidth',1)
title('Normalised voltage magnitude for 30 seconds')
xlabel('time / second')
axis([1,30,0.99,1.01])
ylabel('normalised voltage')
grid on

figure(2)
histfit(Vstd_mov,1200,'GeneralizedExtremeValue');
title('Probability density of standard deviation')
xlabel('standard deviation of voltage in the sliding window')
ylabel('density')
legend('probablity distribution', 'generalized extreme value fit')
grid on;
%std_fit=fitdist(Vstd_mov,'GeneralizedExtremeValue');
