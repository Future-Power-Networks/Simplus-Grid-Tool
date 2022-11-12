% This function plots pole map.

% Author(s): Yitong Li

function PlotPoleMap(EigVecHz,FigN)
  
    figure(FigN);
    
    subplot(1,2,1)
    scatter(real(EigVecHz),imag(EigVecHz),'x','LineWidth',1.5); hold on; grid on;
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
    title('Global pole map');
    
	subplot(1,2,2)
    scatter(real(EigVecHz),imag(EigVecHz),'x','LineWidth',1.5); hold on; grid on;
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
    title('Zoomed pole map');
    axis([-80,20,-150,150]);
end