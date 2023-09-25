%%  tfest function test
load("Dome/vd_rt","vd_rt");
load("Dome/id_rt","id_rt");
load("Dome/vq_rt","vq_rt");
load("Dome/iq_rt","iq_rt");

%u=[vd_rt.vd,vq_rt.vq];
%y=[id_rt.id,iq_rt.iq];
%ymat=tfest(u,y,20,'Ts',1/2400);
P1=timerange(seconds(3),seconds(5-1/fs_new));
P2=timerange(seconds(6),seconds(8-1/fs_new));

vd_rt1=vd_rt(P1,:);
vd_rt2=vd_rt(P2,:);
vq_rt1=vq_rt(P1,:);
vq_rt2=vq_rt(P2,:);
id_rt1=id_rt(P1,:);
id_rt2=id_rt(P2,:);
iq_rt1=iq_rt(P1,:);
iq_rt2=iq_rt(P2,:);

tt=vd_rt1;
tt.vq=vq_rt1.vq;
tt.id=id_rt1.id;
tt.iq=iq_rt1.iq;

sys = tfest(tt,20,'InputName',["vd" "vq"],'OutputName',["id" "iq"]);

%u=[vd_rt.vd,vq_rt.vq];
%y=[id_rt.id,iq_rt.iq];
%ymat=tfest(u,y,20,'Ts',1/2400);

figure(1)
h=bodeplot(sys(1,1)); 
setoptions(h,'FreqUnits','Hz', 'XLim',{[1e-1,1e3]});

% tt=vd_rt;
% tt.id=id_rt.id;
% sys = tfest(tt,20);
% figure(3)
% h=bodeplot(sys); 
% setoptions(h,'FreqUnits','Hz', 'XLim',{[1e-1,1e3]});