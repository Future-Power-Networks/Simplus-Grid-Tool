
%% configuration
%clear;
era_order=12; % order of ERA estimation
f_sample_new=1200; % resample rate
t_start = 3; % start of step in second
t_length =1; % length of data in second
Fs=120000; % original sampling frequency
bus_num=4;
bus_in=2;

%% load data
load("vdq_all_d.mat");
load("vdq_all_q.mat");
load("ListPowerFlowNew.mat");
load("Zsys_SS");

for k=1:bus_num*2
    vdq_d{k}=timeseries2timetable(vdq_all_d{k}.Values);
    vdq_q{k}=timeseries2timetable(vdq_all_q{k}.Values);
end

vdq_d_all=synchronize(vdq_d{1:bus_num*2});
vdq_q_all=synchronize(vdq_q{1:bus_num*2});
vdq_dq_all = synchronize(vdq_d_all, vdq_q_all);


%% trim, resample the data, and remove the steady state
trim = t_start*Fs+1:(t_start+t_length)*Fs;
vdq_dq_trim=retime(vdq_dq_all(trim,:),'regular','SampleRate', f_sample_new);

% removing the steady state value is important for ERA method.
raw_data = vdq_dq_all.Variables;

for i=1:4*bus_num
    vdq_dq_trim(:,i).Variables=vdq_dq_trim(:,i).Variables-raw_data(2.95*Fs,i);
end   
    
figure(111);clf;
plot(vdq_dq_trim.Time,vdq_dq_trim.Data_3_vdq_d_all); title('vdd-bus2');

% ERA and discrete to continuous
clear vdq_dq_s;
for i=1:4*bus_num
    timedata_x=vdq_dq_trim(:,i);
    vdq_dq_s(i,1)=d2c(era(timedata_x, era_order));
end


%% impedance transfer function
v_bus_in = ListPowerFlowNew(bus_in,4);
s=tf('s');
Wbase=2*pi*50;
id_step = -v_bus_in/((13.8626)*10) /s; %-0.02/v_bus2/s;%
iq_step = -v_bus_in * (0.3*0.01) /s; %-0.02/v_bus2/s;%

clear ZsysC2;
for i=1:bus_num
    ZsysC2(i*2-1,1)=vdq_dq_s(i*2-1)/id_step;
    ZsysC2(i*2,1)=vdq_dq_s(i*2)/id_step;
    ZsysC2(i*2-1,2)=vdq_dq_s(i*2-1+8)/iq_step;
    ZsysC2(i*2,2)=vdq_dq_s(i*2+8)/iq_step;
end
%ZsysC2 = minreal(ZsysC2,1e-5);
%% plot;

P=bodeoptions;
P.Grid='on';
P.XLim={[0.1 f_sample_new/2]};
P.FreqUnits='Hz';
P.PhaseWrapping='off';

Bode_enable=1;
if Bode_enable==1
    for bus_k=1:4
        figure(bus_k);clf; 
        subplot(2,2,1);%dd
        bode(ZsysC2(bus_k*2-1,1),P);
        hold on;
        bode(Zsys_SS(bus_k*2-1,bus_in*2-1),P);
        title("Z-dd");
        legend('Estimated from PMU', 'Theorical Results')
        
        subplot(2,2,2); %dq
        bode(ZsysC2(bus_k*2-1,2),P);
        hold on;title("bode-dq");
        bode(Zsys_SS(bus_k*2-1,bus_in*2),P);
        title("Z-dq");
        
        subplot(2,2,3); %qd
        bode(ZsysC2(bus_k*2,1),P);
        hold on;
        bode(Zsys_SS(bus_k*2,bus_in*2-1),P);
        title("Z-qd");
        
        subplot(2,2,4); %qq
        title("bode-qq");
        bode(ZsysC2(bus_k*2,2),P);
        hold on;
        bode(Zsys_SS(bus_k*2,bus_in*2),P);
        title("Z-qq");
    end
end

%run Wholesys_restore.m

