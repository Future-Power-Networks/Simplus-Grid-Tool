t_start = 3; % start time
t_end =5; % end time
fs_new = 1200; % new sample rate

trim = (t_start*Fs+1):t_end*Fs;

vd=out.e_data2{1}.Values; vd0=vd.Data(t_start*Fs-1);
vq=out.e_data2{2}.Values; vq0=vq.Data(t_start*Fs-1);
id=out.e_data2{3}.Values; id0=id.Data(t_start*Fs-1);
iq=out.e_data2{4}.Values; iq0=iq.Data(t_start*Fs-1);

vd_t=timeseries2timetable(vd);
vq_t=timeseries2timetable(vq);
id_t=timeseries2timetable(id);
iq_t=timeseries2timetable(iq);

vd_t=vd_t(trim,:); vd_t.vd = vd_t.vd-vd0; % remove steady state
vq_t=vq_t(trim,:); vq_t.vq = vq_t.vq-vq0;
id_t=id_t(trim,:); id_t.id = id_t.id-id0;
iq_t=iq_t(trim,:); iq_t.iq = iq_t.iq-iq0;
vd_rt=retime(vd_t,'regular','SampleRate', fs_new); % reconfigure sample rate
vq_rt=retime(vq_t,'regular','SampleRate', fs_new);
id_rt=retime(id_t,'regular','SampleRate', fs_new);
iq_rt=retime(iq_t,'regular','SampleRate', fs_new);
save("Dome/vd_rt","vd_rt"); % save the data vd, vq, id, iq
save("Dome/vq_rt","vq_rt");
save("Dome/id_rt","id_rt");
save("Dome/iq_rt","iq_rt");

% Gvd=era(vd_rt,25);
% Gid=era(id_rt,25);
% 
% save("Dome/Gvd","Gvd");
% save("Dome/Gid","Gid");
% 
% Gvd_sys=ss(Gvd.A,Gvd.B, Gvd.C, Gvd.D);
% Gvdtf=tf(Gvd_sys);
% 
%  %id_tf=tf(Gid);
% Gidtf=idtf(Gid);
% 
% %Yddtf=Gidtf/Gvdtf;
% %ydd=-Gid/Gvd;
% figure(1); clf; plot(vd_rt.Time, vd_rt.vd)
% figure(2); clf; 
% h=bodeplot(Gvd); 
% setoptions(h,'FreqUnits','Hz', 'XLim',{[1e-1,1e3]});
% hold on;
% h2=bodeplot(Gid);
% hold on;
% h3=bodeplot(Gvdtf);

%figure(3);clf;
%h=bodeplot(Yddtf); 
%setoptions(h,'FreqUnits','Hz', 'XLim',{[1e-1,1e3]});
%bode(Gid/Gvd);

% trim = t_start*Fs+1:(t_start+t_length)*Fs;
% 
% y_data_all = timeseries2timetable(out.y_rec{1}.Values);
% 
% %y_data_trim = y_data_all; %y_data_all(trim,:);
% f_sample_new=1200;
% y_data_rt = retime(y_data_trim,'regular','SampleRate', f_sample_new);
% 
% figure(1);
% plot(y_data_rt.Time,y_data_rt.Data);
% G_era=era(y_data_rt,25);
% figure(2);
% bode(G_era);
% 
% 
% e = eig(G_era.A);
% %%
% figure(3);
% %subplot(1,2,1)
% scatter(real(pole_sys),imag(pole_sys),'x','LineWidth',1.5); hold on; grid on;
% scatter(real(e),imag(e),'o','LineWidth',1.5);
% xlabel('Real Part (Hz)');
% ylabel('Imaginary Part (Hz)');
% title('Global pole map');
% 
% 
% %y_data_trim = y_data_all(t_start*Fs+1:(t_start+1)*Fs);
% %figure(2);
% %plot(y_data_trim); % plot the impulse response
% 
% %G_era = era(y_data_trim);
% 
% 
% 
% % V2 = timeseries2timetable(out.V_log{1}.Values);
% % V2 = V2(5988000:6036000,:);%6120000
% % V2.Properties.VariableNames(1) = "data";
% % V2.data = 1000000*(V2.data - V2.data(1));
% % V2 = V2(12000:end,:);
% % 
% % figure
% % plot(V2.Time,V2.data)
% % hold on
% % 
% % f_sample_new = 1600;
% % V2_rt = retime(V2,'regular','SampleRate',f_sample_new);
% % V2_rt.data(1) = 0;
% % % plot(V2_rt.Time,V2_rt.data)