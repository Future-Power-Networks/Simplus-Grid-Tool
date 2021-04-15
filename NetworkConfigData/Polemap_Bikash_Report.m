    
    figure(1)
    clf
    scatter(real(pole_sys*2*pi),imag(pole_sys*2*pi),'x','LineWidth',1.5); hold on; grid on;
    xlabel('Real Part (rad)');
    ylabel('Imaginary Part (rad)');
    title('Zoomed pole map');
    axis([-3.2,1,-15,15]);