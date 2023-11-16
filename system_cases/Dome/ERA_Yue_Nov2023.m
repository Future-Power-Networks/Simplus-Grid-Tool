clear vd2_rt;
era_order=12;
f_sample_new=1200;

t_start = 3; % start of step
t_length =0.2;

trim = t_start*Fs+1:(t_start+t_length)*Fs;

vd2_all = timeseries2timetable(out.vdq{1}.Values);
vq2_all = timeseries2timetable(out.vdq{2}.Values);

vd2_trim = vd2_all(trim,:);
vq2_trim = vq2_all(trim,:);
vd2_rt = retime(vd2_trim,'regular','SampleRate', f_sample_new);
vq2_rt = retime(vq2_trim,'regular','SampleRate', f_sample_new);

vd2_rt.vd2=vd2_rt.vd2-vd2_all.vd2(2.8*Fs); % remove steady-state: very important for ERA!!
vq2_rt.Data = vq2_rt.Data-vq2_all.Data(2.8*Fs);

%u_amp=0.015*ListPowerFlowNew(2,4);
v_bus2=ListPowerFlowNew(2,4);
u_amp = v_bus2/((13.8626)*Zbase*5); %current step
s = tf('s');
u_tf_step=-u_amp/s;

figure(1);clf;
subplot(2,1,1);
plot(vd2_rt.Time,vd2_rt.vd2);
subplot(2,1,2);
plot(vq2_rt.Time, vq2_rt.Data);

G_vdr=era(vd2_rt,era_order);
G_vqr=era(vq2_rt,era_order);

vd_tf_d=(tf(G_vdr));
vd_tf = d2c(vd_tf_d);
vq_tf_d=tf(G_vqr);
vq_tf=d2c(vq_tf_d);
%u_tf_d=c2d(u_tf,1/f_sample_new,'zoh');

%zdd1=y_tf_d/u_tf_d;%y_tf_d/u_tf_d;
zdd1=vd_tf/u_tf_step;
zqd1=vq_tf/u_tf_step;
%zdd1=y_tf/u_tf_impulse;
%ydd1=u_tf_d/y_tf_d;

 P=bodeoptions;
 P.Grid='on';
 P.XLim={[0.1 f_sample_new/2]};
 P.FreqUnits='Hz';
% figure(2);clf;
% bode(ydd1,P);
% hold on;
% bode(GminSS(6,5),P);

Zsys_SS = SimplusGT.WholeSysZ_cal(GmObj,YbusObj,Port_i,Port_v);
figure(3);clf;
subplot(2,1,1);
bode(zdd1,P);
hold on;
bode(Zsys_SS(3,3),P);
subplot(2,1,2);
bode(zqd1,P);
hold on;
bode(Zsys_SS(4,3));
%bode(GminSS(Port_i(3),Port_v(3)),P);
