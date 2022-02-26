load('Scope_S13_detuned.mat');
load('Scope_S13_tuned.mat');

Ts=4e-5;
t_start = 34.5;
t_end = 36;
d_start = t_start/Ts+1;
d_end = t_end/Ts+1;
tout = Scope_S13_detuned.time;
A13_detuned = Scope_S13_detuned.signals.values;
A13_tuned = Scope_S13_tuned.signals.values;

color1 = [1,0,0];
color2 = [0,0.45,0.74];
color3 = [0,0.5,0];
L_width=1.2;
Pwidth = 0.25*3.5;
Pheight = 0.8;
PInterval = 0.033;
Fsize = 14;
P1 = [0.1, 0.1,    Pwidth,   Pheight];
P2 = [0.1, 0.1+(Pheight)*2+0.05,  Pwidth,   Pheight];
%P3 = [0.1+(Pwidth+PInterval)*2, 0.1,   Pwidth,   Pheight];
width_L=0.5;%2; %change wdith of the picture.

figure(1);
set(gcf,'position',[500,500,900*1.4*width_L,300*1.4]);
%subplot('Position',P1);
subplot(2,1,1)
plot(tout(d_start:d_end), A13_detuned(d_start:d_end),'LineWidth',L_width, 'Color',color1);
axis([t_start,t_end,0,0.25]) ;
grid on;
subplot(2,1,2)
plot(tout(d_start:d_end), A13_tuned(d_start:d_end),'LineWidth',L_width, 'Color',color2);
axis([t_start,t_end,0,0.25]);
%set(gca,'FontSize',Fsize, 'FontName','TimesNewRoman','FontWeight','Bold');
grid on;