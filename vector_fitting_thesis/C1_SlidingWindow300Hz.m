%% Basic info.
clear
fs=20e3; %20 kHz
Ts=1/fs;
Ttotal=30;
N=fs*Ttotal;
t=(0:1/fs:(Ttotal-1/fs));
load Noise_M;
load Idn;
Idn_ac=Idn-mean(Idn);
noise = Idn_ac;%Noise_M;%Idn_ac;%Noise_M;%Idn_ac;%sin(2*pi*300*t);

f_mea=300;
n_cycle=5;
nc=round((fs/f_mea));
Nc=round(n_cycle*(fs/f_mea));

Ai=0.7;

%% Noise sliding window
1i;
Aid_vect=zeros(N-Nc,2);
for j=1:(N-Nc)
    Aid_c(j)=0;
    Aid_s(j)=0;
    Aid_complex(j)=0;
    noise_s=noise(j:j+Nc-1);
    for k=1:Nc % intergration
        Aid_c(j)=Aid_c(j)+(Ai*cos(k*Ts*2*pi*f_mea)+noise_s(k))*cos(k*Ts*2*pi*f_mea)*Ts/(Nc*Ts/2);
        Aid_s(j)=Aid_s(j)+(Ai*cos(k*Ts*2*pi*f_mea)+noise_s(k))*sin(k*Ts*2*pi*f_mea)*Ts/(Nc*Ts/2);
       % Aid_vect(j,:)=[Aid_c(j),Aid_s(j)];
        Aid_complex(j)=Aid_c(j)+1i*Aid_s(j);
        Aid_abs(j)=abs(Aid_complex(j));
        Aid_err(j)= sqrt((Aid_c(j)-Ai)^2+(Aid_s(j)-0)^2)/Ai;
        %Z_id_abs(i)=sqrt(Z_id_cos(i)^2+Z_id_sin(i)^2);
        %if Z_id_sin(i)>0
        %    Z_id_ang(i)= 360*acos(Z_id_cos(i)/Z_id_abs(i))/(2*pi);
        %else
        %    Z_id_ang(i)= -360*acos(Z_id_cos(i)/Z_id_abs(i))/(2*pi);
        %end
    end
end
figure(1);
clf;
subplot(2,1,1);
histogram(Aid_abs,100);
subplot(2,1,2);
histogram(Aid_err,100);

