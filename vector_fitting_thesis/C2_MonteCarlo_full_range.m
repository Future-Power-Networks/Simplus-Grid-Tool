%% Basic info.
clear
tic
fs=20e3; %20 kHz
Ts=1/fs;
n_cycle=5;
load Fmodelling % fm, with each single frequency point 
load sigma2 % power of the white Guassian Noise.
%load('..\2022\vfit3\myfile_0.mat') % the injection data
load('GROUP1. 15V.mat')
f_inj=Bode_0(2,1:2:end);

%% Monte Carlo Method
Monte_N=2000;  % total times for Monte Carlo sampling
ANreal=zeros(Monte_N,length(f_inj));
ANimag=zeros(Monte_N,length(f_inj));
AN_complex=zeros(Monte_N,length(f_inj));
AN_abs=zeros(Monte_N,length(f_inj));
edge_max=0;
for k=1:Monte_N
    for j=1:length(f_inj) % each frequency point
        f_mea=f_inj(j);
        Nc=round(n_cycle*fs/f_mea); % total sampling point for this frequency
        t=0:Ts:(Nc-1)*Ts;
        Noise_M=sqrt(sigma2)*randn(1,Nc); % the random white noise
        theta_i=rand(1,length(fm))*2*pi; % the random angle for each harmonic
        for h=1:length(fm)  % generate the random noise
            Noise_M=Noise_M+fm(h).Id_amp * sin( 2*pi*fm(h).freq_Hz*t+theta_i(h) );
        end
        for p=1:Nc % intergration
            ANreal(k,j)=ANreal(k,j)+(Noise_M(p))*cos(p*Ts*2*pi*f_mea)*Ts/(Nc*Ts/2); % noise after intergration - cos
            ANimag(k,j)=ANimag(k,j)+(Noise_M(p))*sin(p*Ts*2*pi*f_mea)*Ts/(Nc*Ts/2); % noise after intergration - sin
            AN_complex(k,j)=ANreal(k,j)+1i*ANimag(k,j);
            AN_abs(k,j)=abs(AN_complex(k,j));
            if AN_abs(k,j)>edge_max
                edge_max=AN_abs(k,j);
            end
        end
    end
end
toc
%% Histogram analysis
figure(33); %3d histogram.
clf
edges= 0:0.05/100:0.09 ;
Hist_all = histc(AN_abs', edges, 2);
bar3(edges, Hist_all.');

% kernel fit
for j=1:length(f_inj)
    x=AN_abs(:,j);
    pd(j) = fitdist(x,'Kernel');
    ICDF95(j)=pd(j).icdf(0.95);
end

figure(55)
clf
subplot(2,1,1);
x_pdf=0:0.0001:0.1;
z=zeros(length(f_inj),length(x_pdf));
for j=1:1:length(f_inj)
    z(j,:) = pdf(pd(j),x_pdf);
    line(x_pdf,z(j,:))
    hold on
end
title('kernel fit all')
subplot(2,1,2);
x_pdf=0:0.0001:0.1;
z300 = pdf(pd(46),x_pdf); %300Hz
line(x_pdf,z300)
title('kernel fit 300Hz')

figure(66)
clf
line(f_inj,ICDF95) % for 95% probability
title('95% probability of the maximum absolute error')
ylabel('absolute error / A')
xlabel('frequency / Hz')
grid on;

