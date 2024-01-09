%% Calculation of Zsys
% identity column
I_k=zeros(bus_num*2,2); I_k(bus_in*2-1,1)=1; I_k(bus_in*2,2)=1;

%load("YN.mat");
YN=YbusDSS;
%x1=I_k-YN*Zsys_SS(:,3:4);
x1=I_k-YN*ZsysC2;


for i=1:bus_num
    %YAr(i*2-1:i*2,i*2-1:i*2)= x1(2*i-1:2*i,1:2) * inv(Zsys_SS(2*i-1:2*i,3:4));
    YAr(i*2-1:i*2,i*2-1:i*2)= x1(2*i-1:2*i,1:2) / (ZsysC2(2*i-1:2*i,:));
end

ZAr=SimplusGT.DssSwitchInOut(YAr,8);
[~,ZbusDss]=ZbusObj.GetDSS(ZbusObj);
ZbusSS=SimplusGT.dss2ss(ZbusDss);
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
Zpole=pole(Zsys_restore)/2/pi;
%YAr=minreal(YAr,1e-2);

%[~,GmDSS]=GmObj.GetDSS(GmObj);
% [~,YA_full]=GmObj.GetDSS(GmObj);
% YA = YA_full(Port_i,Port_v);

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


