function Layer1PieChart(axes,data,legend)
    p=pie(axes, data);
    pText = findobj(p,'Type','text');
    percentValues = get(pText,'String'); 
    for i=1:length(data)
        if str2double(erase(percentValues{i},'%')) >0.1 
            pText(i).String = strcat(legend(i),' (',percentValues(i), ')');
        else   %for value smaller than 0.1%, don't show a string.
            pText(i).String = '';
        end
    end
    cmap=lines(length(data));
    axes.Colormap=cmap;
end