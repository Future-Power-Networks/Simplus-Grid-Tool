%% Calculation of Zsys
% identity column

I_k=zeros(bus_num*2,2); I_k(bus_in*2-1,1)=1; I_k(bus_in*2,2)=1;

%load("YN.mat");
YN=YbusDSS;
x1_=I_k-YN*Zsys_SS(:,bus_in*2-1:bus_in*2);
x1=I_k-YN*ZsysC2;
bode_cpr_dq(x1_(3:4,1:2),x1(3:4,1:2),50);


for i=1:bus_num
    %YAr(i*2-1:i*2,i*2-1:i*2)= x1(2*i-1:2*i,1:2) * inv(Zsys_SS(2*i-1:2*i,3:4));
    YAr(i*2-1:i*2,i*2-1:i*2)= x1(2*i-1:2*i,1:2) / (ZsysC2(2*i-1:2*i,:));
end
%ZAr=SimplusGT.DssSwitchInOut(YAr,8);
%[~,ZbusDss]=ZbusObj.GetDSS(ZbusObj);
%ZbusSS=SimplusGT.dss2ss(ZbusDss);
%Zsys_restore = feedback(YAr,ZbusSS);
Zsys_restore_ = inv(YN+YAr);%feedback(ZAr,YN);

%lpf_=tf([1],[1/2000,1]);
%Zsys_restore_lpf=Zsys_restore_*lpf_;
%Zsys_restore_min=minreal(Zsys_restore_lpf,1e-3);

%% order reduction
%Zsys_restore=Zsys_restore_min;
Zsys_restore=balred(Zsys_restore_,100);
%Zsys_restore = pade(Zsys_restore_, 20);
%Zsys_restore=minreal(Zsys_restore_lpf);
%Zsys_restore = inv(YN+YAr);
Zpole_Hz=pole(Zsys_restore)/2/pi;
%YAr=minreal(YAr,1e-2);

[~,GmDSS]=GmObj.GetDSS(GmObj);
[~,YA_full]=GmObj.GetDSS(GmObj);
YA = YA_full(Port_i,Port_v);
%bode_cpr_dq(YAr(3:4,3:4),YA(3:4,3:4),60);

%bode_cpr_dq(YAr(3:4,1:2),YA(3:4,1:2),50);

%% pole map
figure(10); clf;
scatter(real(pole_sys),imag(pole_sys),'x','LineWidth',1.5); hold on; grid on;
scatter(real(Zpole_Hz),imag(Zpole_Hz),'o','LineWidth',1.5);
xlabel('Real Part (Hz)');
ylabel('Imaginary Part (Hz)');
title('Zoomed pole map');
axis([-80,20,-150,150]);
plot([-80,0], [-80,0]*10, '--k','LineWidth',1,'Color','blue')
plot([-80,0], [80,0]*10, '--k','LineWidth',1,'Color','blue')
legend('theorical','realisation','10% damping line')

Bode_enable=1;
if Bode_enable==1
    for bus_k=1:4
        figure(bus_k+20);clf; 
        subplot(2,2,1);%dd
        bode(Zsys_restore(bus_k*2-1,bus_k*2-1),P);
        hold on;
        bode(Zsys_SS(bus_k*2-1,bus_k*2-1),P);
        title("Z-dd");
        legend('Estimated', 'Theorical')
        
        subplot(2,2,2); %dq
        bode(Zsys_restore(bus_k*2-1,bus_k*2),P);
        hold on;title("bode-dq");
        bode(Zsys_SS(bus_k*2-1,bus_k*2),P);
        title("Z-dq");
        
        subplot(2,2,3); %qd
        bode(Zsys_restore(bus_k*2,bus_k*2-1),P);
        hold on;
        bode(Zsys_SS(bus_k*2,bus_k*2-1),P);
        title("Z-qd");
        
        subplot(2,2,4); %qq
        title("bode-qq");
        bode(Zsys_restore(bus_k*2,bus_k*2),P);
        hold on;
        bode(Zsys_SS(bus_k*2,bus_k*2),P);
        title("Z-qq");
    end
end

function bode_cpr_dq(z1,z2,figk)
    P=bodeoptions;
    P.Grid='on';
    P.XLim={[0.1 600]};
    P.FreqUnits='Hz';
    P.PhaseWrapping='off';

    figure(figk);clf; 
    subplot(2,2,1);%dd
    bode(z1(1,1),P);
    hold on;
    bode(z2(1,1),P);
    title("Z-dd");
    legend('Estimated', 'Theorical')
    
    subplot(2,2,2); %dq
    bode(z1(1,2),P);
    hold on;title("bode-dq");
    bode(z2(1,2),P);
    title("Z-dq");
    
    subplot(2,2,3); %qd
    bode(z1(2,1),P);
    hold on;
    bode(z2(2,1),P);
    title("Z-qd");
    
    subplot(2,2,4); %qq
    title("bode-qq");
    bode(z1(2,2),P);
    hold on;
    bode(z2(2,2),P);
    title("Z-qq");
 end
