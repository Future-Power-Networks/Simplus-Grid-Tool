clear y_data_rt;

t_start = 3; % start of impluse
t_length =1;

trim = t_start*Fs+1:(t_start+t_length)*Fs;

y_data_all = timeseries2timetable(out.vd_rec{1}.Values);

y_data_trim = y_data_all(trim,:);
f_sample_new=1200;
y_data_rt = retime(y_data_trim,'regular','SampleRate', f_sample_new);

u_amp=1/((13.8626)*Zbase*5); %current impulse amplitude: 1/(10*R）
s = tf('s');
u_tf=u_amp*(1-exp(1)^(-0.1)); % 0.1s duration

%u_tf=1*u_amp/s;

figure(1);clf;
plot(y_data_rt.Time,y_data_rt.Data);


G_era=era(y_data_rt,20);
y_tf_d=(tf(G_era));
y_tf = d2c(y_tf_d);
%u_tf_d=c2d(u_tf,1/f_sample_new,'zoh');

zdd1=y_tf/u_tf;%y_tf_d/u_tf_d;

P=bodeoptions;
P.Grid='on';
P.XLim={[0.1 1e3]};
P.FreqUnits='Hz';
figure(2);clf;
bode(zdd1,P);

Zsys_SS = SimplusGT.WholeSysZ_cal(GmObj,YbusObj,Port_i,Port_v);
figure(3);clf;
bode(Zsys_SS(3,3),P);

%e = eig(G_era.A);
%%
%figure(3);
%subplot(1,2,1)
%scatter(real(pole_sys),imag(pole_sys),'x','LineWidth',1.5); hold on; grid on;
%scatter(real(e),imag(e),'o','LineWidth',1.5);
%xlabel('Real Part (Hz)');
%ylabel('Imaginary Part (Hz)');
%title('Global pole map');


%y_data_trim = y_data_all(t_start*Fs+1:(t_start+1)*Fs);
%figure(2);
%plot(y_data_trim); % plot the impulse response

%G_era = era(y_data_trim);



% V2 = timeseries2timetable(out.V_log{1}.Values);
% V2 = V2(5988000:6036000,:);%6120000
% V2.Properties.VariableNames(1) = "data";
% V2.data = 1000000*(V2.data - V2.data(1));
% V2 = V2(12000:end,:);
% 
% figure
% plot(V2.Time,V2.data)
% hold on
% 
% f_sample_new = 1600;
% V2_rt = retime(V2,'regular','SampleRate',f_sample_new);
% V2_rt.data(1) = 0;
% % plot(V2_rt.Time,V2_rt.data)