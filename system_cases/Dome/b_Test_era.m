%%  tfest function test
load("Dome/vd_rt","vd_rt");
% load("Dome/id_rt","id_rt");
% load("Dome/vq_rt","vq_rt");
% load("Dome/iq_rt","iq_rt");

%u=[vd_rt.vd,vq_rt.vq];
%y=[id_rt.id,iq_rt.iq];
%ymat=tfest(u,y,20,'Ts',1/2400);
P1=timerange(seconds(3),seconds(5-1/fs_new));
P2=timerange(seconds(6),seconds(8-1/fs_new));

vd_rt1=vd_rt(P1,:);
% vd_rt2=vd_rt(P2,:);
% vq_rt1=vq_rt(P1,:);
% vq_rt2=vq_rt(P2,:);
% id_rt1=id_rt(P1,:);
% id_rt2=id_rt(P2,:);
% iq_rt1=iq_rt(P1,:);
% iq_rt2=iq_rt(P2,:);

%% ERA
era_order=10;
vd1=era(vd_rt1,era_order)*2*13.86*100;
% vq1=era(vq_rt1,era_order);
% id1=era(id_rt1,era_order);
% iq1=era(iq_rt1,era_order);
% 
% vd2=era(vd_rt2,era_order);
% vq2=era(vq_rt2,era_order);
% id2=era(id_rt2,era_order);
% iq2=era(iq_rt2,era_order);

% Vsys=[(vd1),(vd2); (vq1),(vq2)];
% Isys=[(id1),(id2);(iq1),(iq2)];
% Ysys_est=Isys/Vsys;

%% figure
figure(1); clf;
h=bodeplot(vd1(1,1)); 
setoptions(h,'FreqUnits','Hz', 'XLim',{[1e-1,1e3]});
hold on;
h2=bodeplot(GminSS(6,5));

%figure(2);clf;
%plot(vd_rt1.Time,[vd_rt1.vd,vd_rt2.vd]);

% figure(3);
% pole_sys=eig(Ysys_est(1,1).A);
% scatter(real(pole_sys),imag(pole_sys),'x','LineWidth',1.5); hold on; grid on;
% %scatter(real(e),imag(e),'o','LineWidth',1.5);
% xlabel('Real Part (Hz)');
% ylabel('Imaginary Part (Hz)');
% title('Global pole map');

% tt=vd_rt;
% tt.id=id_rt.id;
% sys = tfest(tt,20);
% figure(3)
% h=bodeplot(sys); 
% setoptions(h,'FreqUnits','Hz', 'XLim',{[1e-1,1e3]});