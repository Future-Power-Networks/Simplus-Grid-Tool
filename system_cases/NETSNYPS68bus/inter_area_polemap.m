%pole_sys_tuned_inter = pole_sys;
%save pole_sys_tuned_inter
figure(5);
clf
load pole_sys_tuned_inter.mat;
load pole_sys_detuned.mat;

scatter(real(pole_sys_detuned),imag(pole_sys_detuned),'x','LineWidth',1.5); 
hold on; grid on;
xlabel('Real Part (Hz)');
ylabel('Imaginary Part (Hz)');
title('Zoomed pole map');
axis([-0.2,0.02,-1,1]);

hold on;
scatter(real(pole_sys_tuned_inter),imag(pole_sys_tuned_inter),'x','LineWidth',1.5);