%% Detuned system
clear figure 1

sel = 1; % 1 for detuned, 2 for tuned, 3 for tuned against.
if sel == 1
    load('StepR1.mat','out');
elseif sel == 2
    load('StepR2.mat','out');
elseif sel == 3
    load('StepR3.mat','out');
end;
Ts=1/3e4;
t_start = 9.5;
t_end = 16;
d_start = t_start/Ts+1;
d_end = t_end/Ts+1;
tout = out.tout;
A11_S = out.ScopeData_ApparentPower.signals(1).values;
A15_S = out.ScopeData_ApparentPower.signals(2).values;
A29_S = out.ScopeData_ApparentPower.signals(3).values;
color1 = [1,0,0];
color2 = [0,0.45,0.74];
color3 = [0,0.5,0];
L_width=1.7;
Pwidth = 0.25;
Pheight = 0.8;
PInterval = 0.033;
Fsize = 14;
P1 = [0.1, 0.1,    Pwidth,   Pheight];
P2 = [0.1+Pwidth+PInterval, 0.1,   Pwidth,   Pheight];
P3 = [0.1+(Pwidth+PInterval)*2, 0.1,   Pwidth,   Pheight];
width_L=1;%2; %change wdith of the picture.

figure(1);
set(gcf,'position',[500,500,900*1.4*width_L,200*1.4]);
subplot('Position',P1);
plot(tout(d_start:d_end), A11_S(d_start:d_end),'LineWidth',L_width, 'Color',color1);
if sel==3
    axis([t_start,t_end,10.2,11.6]) ;
else
    axis([t_start,t_end,10.8,11.4]);
end
%title('A11 Apparent Power (pu)')
set(gca,'FontSize',Fsize, 'FontName','TimesNewRoman','FontWeight','Bold');
grid on;

subplot('Position',P2);
plot( tout(d_start:d_end), A15_S(d_start:d_end),'LineWidth',L_width,'Color',color2);
if sel == 3
    axis([t_start,t_end,-10,80]);
else
    axis([t_start,t_end,2.5,3.2]);
end
%title('A15 Apparent Power (pu)')
set(gca,'FontSize',Fsize, 'FontName','TimesNewRoman','FontWeight','Bold');
grid on;

subplot('Position',P3);
plot(tout(d_start:d_end), A29_S(d_start:d_end),'LineWidth',L_width, 'Color',color3);
if sel== 3
    axis([t_start,t_end,3.14,3.28]);
else
    axis([t_start,t_end,3.18,3.28]);
end
%title('A29 Apparent Power (pu)')
set(gca,'FontSize',Fsize, 'FontName','TimesNewRoman','FontWeight','Bold');
grid on;

% %% Tuned following the grey box
% load('StepR2.mat','out');
% d_start = t_start/Ts+1;
% d_end = t_end/Ts+1;
% tout = out.tout;
% A11_S = out.ScopeData_ApparentPower.signals(1).values;
% A15_S = out.ScopeData_ApparentPower.signals(2).values;
% A29_S = out.ScopeData_ApparentPower.signals(3).values;
% 
% L_width=0.2;
% figure(1);
% subplot(3,3,4)
% plot(tout(d_start:d_end), A11_S(d_start:d_end),'LineWidth',L_width, 'Color',color1);
% axis([t_start,t_end,10.8,11.4]);
% grid on;
% 
% subplot(3,3,5)
% plot( tout(d_start:d_end), A15_S(d_start:d_end),'LineWidth',L_width,'Color',color2);
% axis([t_start,t_end,2.5,3.2]);
% grid on;
% 
% subplot(3,3,6)
% plot(tout(d_start:d_end), A29_S(d_start:d_end),'LineWidth',L_width, 'Color',color3);
% axis([t_start,t_end,3.18,3.28]);
% grid on;
% 
% 
% 
% %% Tuned against
% 
% load('StepR3.mat','out');
% d_start = t_start/Ts+1;
% d_end = t_end/Ts+1;
% tout = out.tout;
% A11_S = out.ScopeData_ApparentPower.signals(1).values;
% A15_S = out.ScopeData_ApparentPower.signals(2).values;
% A29_S = out.ScopeData_ApparentPower.signals(3).values;
% 
% L_width=0.2;
% figure(1);
% subplot(3,3,7)
% plot(tout(d_start:d_end), A11_S(d_start:d_end),'LineWidth',L_width, 'Color',color1);
% axis([t_start,t_end,10.2,11.6]) ;
% grid on;
% 
% subplot(3,3,8)
% plot( tout(d_start:d_end), A15_S(d_start:d_end),'LineWidth',L_width,'Color',color2);
% axis([t_start,t_end,-10,80]);
% grid on;
% 
% subplot(3,3,9)
% plot(tout(d_start:d_end), A29_S(d_start:d_end),'LineWidth',L_width, 'Color',color3);
% axis([t_start,t_end,3.14,3.28]);
% grid on;
% 
% 
% %%
% subplot(3,3,1)
% title('A11 Apparent Power (pu)')
% subplot(3,3,2)
% title('A15 Apparent Power (pu)')
% subplot(3,3,3)
% title('A29 Apparent Power (pu)')
% 
% subplot(3,3,4)
% title('A11 Apparent Power (pu)')
% subplot(3,3,5)
% title('A15 Apparent Power (pu)')
% subplot(3,3,6)
% title('A29 Apparent Power (pu)')
% 
% subplot(3,3,7)
% title('A11 Apparent Power (pu)')
% subplot(3,3,8)
% title('A15 Apparent Power (pu)')
% subplot(3,3,9)
% title('A29 Apparent Power (pu)')
