fig = figure(2);
clf;
fig.Position = [600 300 1000 700];
h=scatter(real(pole_sys),imag(pole_sys),'x','LineWidth',1.5); hold on; grid on;
xlabel('Real Part (Hz)');
ylabel('Imaginary Part (Hz)');
title('Zoomed pole map');
axis([-20,2,-100,100]);
plot([-80,0], [-80,0]*10, '--k','LineWidth',2,'Color','blue');
plot([-80,0], [80,0]*10, '--k','LineWidth',2,'Color','blue');
%legend('mode','10% damping line');
h.ButtonDownFcn = @ModeSelect;
%b=uibutton(fig);

function [coordinateSelected, minIdx]= ModeSelect(hObj, event)
x = get(hObj, 'XData');
y = get(hObj, 'YData');
pt = event.IntersectionPoint(1:2);       % The (x0,y0) coordinate just selected
coordinates = [x(:),y(:)];     % matrix of your input coordinates
dist_m = pdist2(pt,coordinates);  
[~, minIdx] = min(dist_m);            % index of minimum distance to points
coordinateSelected = coordinates(minIdx,:); %the selected coordinate
% --- Get index of the clicked point
%[~, i] = min((e.IntersectionPoint(1)-x).^2 + (e.IntersectionPoint(2)-y).^2);
%plot(x(minIdx), y(minIdx), 'o', 'color', 'red');
fprintf('the selected mode is %.3f + j%.3f Hz, its index is %d\n', x(minIdx), y(minIdx), minIdx);
end

