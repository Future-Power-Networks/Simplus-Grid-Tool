function figo = PlotHeatMap(x,y,v,fig, climt)

    % Find the figure
    try
        if(isnumeric(fig))
            fig = figure(fig);
        end
    catch
        fig = figure();
    end
    figo = fig;

    % Interpolant
    F = scatteredInterpolant(x,y,v);
    [xq,yq] = meshgrid(-4.5:0.1:4.5);
    % F.Method = 'linear';
    F.Method = 'natural';
    vq = F(xq,yq);
    
    % Further deal with unreasonable data
%     [rmax,cmax] = size(vq);
%     for r = 1:rmax
%         for c = 1:cmax
%             if vq(r,c) < 0
%                 vq(r,c) = 0;
%             elseif vq(r,c) > max(v)
%                 vq(r,c) = max(v);
%             end
%         end
%     end
    
    % Plot
    figure(fig);
    contourf(xq,yq,vq,150,'LineColor','none');      % Color map
    colorbar;                                       % Color bar
    
    %if ((min(v)<range(1)) || (max(v)>range(2)) || (range(1)>range(2)))
    %    error('range is not set properly.');
    %end
    

    % red - yellow - green - white
    c1 = [200, 0, 0]/255;
    c2 = [255, 192, 0]/255;
    c3 = [255, 255, 89]/255;
    c4 = [147, 255, 196]/255;
    c5 = [213, 249, 243]/255;
    c6 = [255,255,255]/255;
    
    % red - yellow - blue - white
%     c1 = [200, 0, 0]/255;
%     c2 = [255, 192, 0]/255;
%     c3 = [255, 255, 89]/255;
%     c4 = [132, 221, 255]/255;
%     c5 = [182, 235, 255]/255;
%     c6 = [255,255,255]/255;
    
    for i=1:3
         g1=linspace(c1(i),c2(i),200);
         g2=linspace(c2(i),c3(i),300);
         g3=linspace(c3(i),c4(i),200);
         g4=linspace(c4(i),c5(i),200);
         g5=linspace(c5(i),c6(i),500);
        gx(:,i)=[g1,g2,g3,g4,g5]';
    end
    %colormap(gx);
    %fig.Colormap = colormap(flipud(jet));%[GradRed GradGreen GradBlue];
    fig.Colormap = colormap(gx);
    clim(climt)
   	%figure;
    %mesh(xq,yq,vq);
    %hold on;
    %plot3(x,y,v,'.');
    %set(gca,'ColorScale','log')
end

function m = Affine(x,y,r)
    m = x*(1-r) + y*r;
end