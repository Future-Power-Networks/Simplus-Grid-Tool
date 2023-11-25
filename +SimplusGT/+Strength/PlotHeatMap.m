% Plot heat map based on a layout graph and its x, y, z data.
%
% Author(s): Yitong, Yue Zhu

function PlotHeatMap(x,y,z)

    % Interpolant
    Delta = max(max(x)-min(x),max(y)-min(y));
    xyLimit = max(max(abs(x)),max(abs(y))) + Delta*0.1;
    [xmesh,ymesh] = meshgrid(linspace(-xyLimit,xyLimit,100));
    
    F = scatteredInterpolant(x,y,z);
    F.Method = 'natural';
    zmesh = F(xmesh,ymesh);

    % Plot
    contourf(xmesh,ymesh,zmesh,150,'LineColor','none');      % Color map
    
    % Set color bar
    colorbar;                                       % Color bar
    
    % red - yellow - green - white
    c1 = [200, 0, 0]/255;
    c2 = [255, 192, 0]/255;
    c3 = [255, 255, 89]/255;
    c4 = [147, 255, 196]/255;
    c5 = [213, 249, 243]/255;
    c6 = [255,255,255]/255;
    
    for i=1:3
         g1=linspace(c1(i),c2(i),200);
         g2=linspace(c2(i),c3(i),300);
         g3=linspace(c3(i),c4(i),200);
         g4=linspace(c4(i),c5(i),200);
         g5=linspace(c5(i),c6(i),500);
        gx(:,i)=[g1,g2,g3,g4,g5]';
    end

    colormap(gx);

    % Limit the graph z data
    climt = [min(z), max(z)];
    clim(climt);

end