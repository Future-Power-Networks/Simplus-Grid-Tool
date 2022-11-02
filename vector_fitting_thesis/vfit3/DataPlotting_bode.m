load('myfile_0.mat')
p_start=10;
%% 
Length=length(Bode_0);
Frequency=Bode_0(2,1:2:Length-1);

Ydd_mag=Bode_0(3,1:2:Length-1);
Ydd_mag_db = 20*log10(Ydd_mag); 
Ydd_Phase=Bode_0(4,1:2:Length-1);

Ydq_mag=Bode_0(3,2:2:Length);
Ydq_mag_db = 20*log10(Ydq_mag);
Ydq_Phase=Bode_0(4,2:2:Length);

Yqd_mag=Bode_0(5,1:2:Length-1);
Yqd_mag_db = 20*log10(Yqd_mag);
Yqd_Phase=Bode_0(6,1:2:Length-1);

Yqq_mag=Bode_0(5,2:2:Length);
Yqq_mag_db = 20*log10(Yqq_mag);
Yqq_Phase=Bode_0(6,2:2:Length);

Length_new=Length/2;
%% Plot
% Arrange x axis
fbd_L = min(Frequency);
fbd_H = max(Frequency);
xmax = log10(fbd_H)-log10(fbd_L)+1;
x_v(1) = fbd_L;
for x=2:(xmax-1)
x_v(x) = fbd_L*(10^(x-1));
end
x_v(xmax) = fbd_H;

% Color
color_blue = [0, 0.4470, 0.7410];
color2 = [0.8500    0.3250    0.0980];

Mag_H = -20;
Mag_L = -35;
Mag_v = [-35,-30,-25,-20];
Pha_H = 60;
Pha_L = -30;
Pha_v = [-30,-15,0,15,30,45,60];

%% dd
figure(p_start+1)
subplot(2,1,1);
semilogx(Frequency,Ydd_mag_db,'linewidth',1.5,'color',color_blue); grid on; hold on;
semilogx(Frequency,Ydd_mag_db,'rx','linewidth',1,'color',color2,'MarkerSize',8); grid on; hold off;
set(gca,'YLim',[Mag_L Mag_H]);
set(gca,'YTick',Mag_v);
set(gca,'XLim',[fbd_L fbd_H]);
set(gca,'XTick',x_v);
ylabel('Magnitude (dB)')

subplot(2,1,2);
semilogx(Frequency,Ydd_Phase,'linewidth',1.5,'color',color_blue); grid on; hold on;
semilogx(Frequency,Ydd_Phase,'rx','linewidth',1,'color',color2,'MarkerSize',8); grid on; hold off;
set(gca,'YLim',[Pha_L Pha_H]);
set(gca,'YTick',Pha_v);
set(gca,'XLim',[fbd_L fbd_H]);
set(gca,'XTick',x_v);
ylabel('Phase (degree)')
xlabel('Frequency (Hz)')
figure(p_start+1); mtit('$Y_{dd}$','interpreter','latex','FontSize',12)   

%% dq
figure(p_start+2)
subplot(2,1,1);
semilogx(Frequency,Ydq_mag_db,'linewidth',1.5,'color',color_blue); grid on; hold on;
semilogx(Frequency,Ydq_mag_db,'rx','linewidth',1,'color',color2,'MarkerSize',8); grid on; hold off;
set(gca,'YLim',[Mag_L Mag_H]);
set(gca,'YTick',Mag_v);
set(gca,'XLim',[fbd_L fbd_H]);
set(gca,'XTick',x_v);
ylabel('Magnitude (dB)')

subplot(2,1,2);
semilogx(Frequency,Ydq_Phase,'linewidth',1.5,'color',color_blue); grid on; hold on;
semilogx(Frequency,Ydq_Phase,'rx','linewidth',1,'color',color2,'MarkerSize',8); grid on; hold off;
set(gca,'YLim',[Pha_L Pha_H]);
set(gca,'YTick',Pha_v);
set(gca,'XLim',[fbd_L fbd_H]);
set(gca,'XTick',x_v);
ylabel('Phase (degree)')
xlabel('Frequency (Hz)')
figure(p_start+2); mtit('$Y_{dq}$','interpreter','latex','FontSize',12)   

%% qd
figure(p_start+3)
subplot(2,1,1);
semilogx(Frequency,Yqd_mag_db,'linewidth',1.5,'color',color_blue); grid on; hold on;
semilogx(Frequency,Yqd_mag_db,'rx','linewidth',1,'color',color2,'MarkerSize',8); grid on; hold off;
set(gca,'YLim',[Mag_L Mag_H]);
set(gca,'YTick',Mag_v);
set(gca,'XLim',[fbd_L fbd_H]);
set(gca,'XTick',x_v);
ylabel('Magnitude (dB)')

subplot(2,1,2);
semilogx(Frequency,Yqd_Phase,'linewidth',1.5,'color',color_blue); grid on; hold on;
semilogx(Frequency,Yqd_Phase,'rx','linewidth',1,'color',color2,'MarkerSize',8); grid on; hold off;
set(gca,'YLim',[Pha_L Pha_H]);
set(gca,'YTick',Pha_v);
set(gca,'XLim',[fbd_L fbd_H]);
set(gca,'XTick',x_v);
ylabel('Phase (degree)')
xlabel('Frequency (Hz)')
figure(p_start+3); mtit('$Y_{qd}$','interpreter','latex','FontSize',12)   

%% qq
figure(p_start+4)
subplot(2,1,1);
semilogx(Frequency,Yqq_mag_db,'linewidth',1.5,'color',color_blue); grid on; hold on;
semilogx(Frequency,Yqq_mag_db,'rx','linewidth',1,'color',color2,'MarkerSize',8); grid on; hold off;
set(gca,'YLim',[Mag_L Mag_H]);
set(gca,'YTick',Mag_v);
set(gca,'XLim',[fbd_L fbd_H]);
set(gca,'XTick',x_v);
ylabel('Magnitude (dB)')

subplot(2,1,2);
semilogx(Frequency,Yqq_Phase,'linewidth',1.5,'color',color_blue); grid on; hold on;
semilogx(Frequency,Yqq_Phase,'rx','linewidth',1,'color',color2,'MarkerSize',8); grid on; hold off;
set(gca,'YLim',[Pha_L Pha_H]);
set(gca,'YTick',Pha_v);
set(gca,'XLim',[fbd_L fbd_H]);
set(gca,'XTick',x_v);
ylabel('Phase (degree)')
xlabel('Frequency (Hz)')
figure(p_start+4); mtit('$Y_{qq}$','interpreter','latex','FontSize',12)   
