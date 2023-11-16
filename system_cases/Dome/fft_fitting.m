clear y_data_rt;

t_start = 3; % start of impluse
t_length =1;

trim = t_start*Fs+1:(t_start+t_length)*Fs;
y_data_all = timeseries2timetable(out.vd2{1}.Values);
y_data_trim = y_data_all(trim,:);
f_sample_new=1200;
y_data_rt = retime(y_data_trim,'regular','SampleRate', f_sample_new);

figure(1);clf;
plot(y_data_rt.Time,y_data_rt.Data);

Y=fft(y_data_rt.Data);

fx = Fs*(0:(n/2))/n;
P = abs(Y/n).^2;