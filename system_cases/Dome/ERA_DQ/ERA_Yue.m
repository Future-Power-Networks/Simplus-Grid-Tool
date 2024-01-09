
%% configuration
era_order=12; % order of ERA estimation
f_sample_new=1200; % resample rate
t_start = 3; % start of step in second
t_length =1; % length of data in second
Fs=120000; % original sampling frequency
%% load data
load("vdq2_d.mat");
load("vdq2_q.mat");
load("vdq3_d.mat");
load("vdq3_q.mat");
load("ListPowerFlowNew.mat");

for k=2:3
    if k==2
        vd_d=timeseries2timetable(vdq2_d{1}.Values);
        vq_d=timeseries2timetable(vdq2_d{2}.Values);
        vd_q=timeseries2timetable(vdq2_q{1}.Values);
        vq_q=timeseries2timetable(vdq2_q{2}.Values);
    elseif k==3
        vd_d=timeseries2timetable(vdq3_d{1}.Values);
        vq_d=timeseries2timetable(vdq3_d{2}.Values);
        vd_q=timeseries2timetable(vdq3_q{1}.Values);
        vq_q=timeseries2timetable(vdq3_q{2}.Values);
    end
    
    %% trim, resample the data, and remove the steady state
    trim = t_start*Fs+1:(t_start+t_length)*Fs;
    vdd_trim = retime(vd_d(trim,:),'regular','SampleRate', f_sample_new);
    vqd_trim = retime(vq_d(trim,:),'regular','SampleRate', f_sample_new);
    vdq_trim = retime(vd_q(trim,:),'regular','SampleRate', f_sample_new);
    vqq_trim = retime(vq_q(trim,:),'regular','SampleRate', f_sample_new);
    % removing the steady state value is important for ERA method.
    vdd_trim.Data=vdd_trim.Data-vd_d.Data(2.8*Fs);
    vqd_trim.Data=vqd_trim.Data-vq_d.Data(2.8*Fs);
    vdq_trim.Data=vdq_trim.Data-vd_q.Data(2.8*Fs);
    vqq_trim.Data=vqq_trim.Data-vq_q.Data(2.8*Fs);
    
    
    figure(k*10);clf;
    title("Data for ERA")
    subplot(2,2,1);
    plot(vdd_trim.Time,vdd_trim.Data); title('vdd');
    subplot(2,2,2);
    plot(vdq_trim.Time,vdq_trim.Data);title('vdq');
    subplot(2,2,3);
    plot(vqd_trim.Time,vqd_trim.Data);title('vqd');
    subplot(2,2,4);
    plot(vqq_trim.Time,vqq_trim.Data);title('vqq');
    
    %% ERA estimation, discrete to continuous
    vdd_s=d2c(era(vdd_trim, era_order));
    vqd_s=d2c(era(vqd_trim, era_order));
    vdq_s=d2c(era(vdq_trim, era_order));
    vqq_s=d2c(era(vqq_trim, era_order));
    
    %% impedance transfer function
    v_bus2 = ListPowerFlowNew(2,4);
    s=tf('s');
    Wbase=2*pi*50;
    id_step = -v_bus2/((13.8626)*10) /s; %-0.02/v_bus2/s;%
    iq_step = -v_bus2 * (0.3*0.01) /s; %-0.02/v_bus2/s;%
    
    zdd=vdd_s/id_step;
    zqd=vqd_s/id_step;
    zdq=vdq_s/iq_step;
    zqq=vqq_s/iq_step;

    Zsys_r{k,2}=[zdd,zdq;zqd,zqq];
        
    
    %% plot;
    load("Zsys_SS");
    P=bodeoptions;
    P.Grid='on';
    P.XLim={[0.1 f_sample_new/2]};
    P.FreqUnits='Hz';
    
    figure(k);clf; 
    subplot(2,2,1);%dd
    bode(zdd,P);
    hold on;
    bode(Zsys_SS(3,k*2-1),P);
    title("Z-dd");
    legend('Estimated from PMU', 'Theorical Results')
    
    subplot(2,2,2); %dq
    bode(zdq,P);
    hold on;title("bode-dq");
    bode(Zsys_SS(3,k*2),P);
    title("Z-dq");
    
    subplot(2,2,3); %qd
    bode(zqd,P);
    hold on;
    bode(Zsys_SS(4,k*2-1),P);
    title("Z-qd");
    
    subplot(2,2,4); %qq
    title("bode-qq");
    bode(zqq,P);
    hold on;
    bode(Zsys_SS(4,k*2),P);
    title("Z-qq");
end
