function Layer2BarChart(ax,data,legend)
    
    X = categorical(legend);
    X = reordercats(X,legend);
    b=bar(ax,X, data);
    
    bar_num=length(data);
    cmap=lines(bar_num);
    b.FaceColor='flat';
    for i=1:bar_num
        b.CData(i,:)=cmap(i,:);
    end
%    b.CData(1,:) = [1 0 0];
%    b.CData(2,:) = [0 1 0];
%    b.CData(3,:) = [0 0 1];
end