clear
Ydmag=[];
Ydang=[];
for i=6:15   % file 6 to 15.
    file_name=strcat('Injection300Hz_7Vamp/myfile_0',num2str(i),'.mat'); %inject amplitude: 15V, not 7V, 7V is just the folder name.
    load(file_name)
    Yd_mag_new=Bode_0(3,1:2:length(Bode_0));
    Yd_angle_new=(Bode_0(4,1:2:length(Bode_0)))/360*2*pi; % correct phase: -90 degree!
    Ydmag=[Ydmag,Yd_mag_new];
    Ydang=[Ydang,Yd_angle_new];
end

figure(99)
histogram(Ydmag,40);
for k=1:length(Ydmag)
    Yd_complex(k)=Ydmag(k)*exp(1i*Ydang(k));
end
%Yd_complex=Ydmag*exp(1i*Ydang);
figure(88)
clf
Yd_real=real(Yd_complex);
Yd_imag=imag(Yd_complex);
plot(Yd_real);
hold on;
plot(Yd_imag);

%Yd_Mag_true=mean(Ydmag);
Yd_real_true=mean(Yd_real);
Yd_imag_true=mean(Yd_imag);


for i=1:length(Ydmag)
    %error(i)=abs(Ydmag(i)-Yd_Mag_true)/Yd_Mag_true;
    error(i)= sqrt((Yd_real(i)-Yd_real_true)^2+(Yd_imag(i)-Yd_imag_true)^2)/...
        sqrt(Yd_real_true^2+Yd_imag_true^2);
end
max(error)

figure(100)
clf;
histogram(error,80);
grid on;
title('Histogram of relative errors from 1000 tests')
ylabel('appeared times')
xlabel('relative error')