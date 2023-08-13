
%% Load map data
% DataName = 'K_68Bus_IBR_Load_Data';
% DataName = 'K_68Bus_IBR_Data';
% DataName = 'K_68Bus_IBR_17_Data';
% DataName = 'K_68Bus_IBR_17_14_Data';
DataName = 'K_68Bus_IBR_17_14_7_Data';

Data = load(DataName).SaveData;

%%
KH              = Data.KH;
YbusVI          = Data.YbusVI;
YbusVIF         = Data.YbusVIF;
GbusVI          = Data.GbusVI;
GbusVIF         = Data.GbusVIF;
YbusOrigin      = Data.YbusOrigin;
%Index_Vbus      = Data.Index_Vbus;
%Index_Ibus      = Data.Index_Ibus;
%Index_Fbus      = Data.Index_Fbus;
%Index_Ebus      = Data.Index_Ebus;
Order_Old2New   = Data.Order_Old2New;
Order_New2Old   = Data.Order_New2Old;

%%
FigNum = 0;
ColorRGB();
FigSize = [0.1 0.1 0.5 0.75];

FigNum = FigNum + 1;
figure(FigNum)
set(gcf,'units','normalized','outerposition',FigSize);

%% Plot graph
GraphMatrix = NormMatrixElement(YbusOrigin,'DiagFlag',0);
GraphData = graph(GraphMatrix,'upper');
GraphFigure = plot(GraphData); grid on; hold on;
highlight(GraphFigure,GraphData,'EdgeColor',[0,0,0],'LineWidth',1.1);       % Change all edges to black by default
highlight(GraphFigure,GraphData,'NodeColor',[0,0,0]);                    	% Change all nodes to black by default
highlight(GraphFigure,GraphData,'MarkerSize',9);
highlight(GraphFigure,GraphData,'NodeFontSize',12);
highlight(GraphFigure,GraphData,'NodeFontWeight','bold');

%% sort out SG-bus, IBR-bus and floating bus
k1=1;
k2=1;
k3=1;
for k=1:N_Bus
    if ApparatusType{k} >= 10 && ApparatusType{k} <20 % GFL
        Index_Ibus(k1)=k;
        k1=k1+1;
    elseif (ApparatusType{k} >= 20 && ApparatusType{k} <40) || ApparatusType{k} <10% GFM & SG
        Index_Vbus(k2)=k;
        k2=k2+1;
    elseif ApparatusType{k}==100
        Index_Fbus(k3)=k;
        k3=k3+1;
    end
end

%% Set current node
highlight(GraphFigure,Index_Ibus,'NodeColor',[0,128,0]/255);          % Change all current node to green by default

%% Set voltage node
highlight(GraphFigure,Index_Vbus,'NodeColor',[0,0,0]);      	% Change all voltage node to black by default

%% Set floating node
highlight(GraphFigure,Index_Fbus,'NodeColor',[0.7,0.7,0.7]);   	% Change all floating node to gray by default


XData_ = GraphFigure.XData';
YData_ = GraphFigure.YData';

XData = XData_([CIMR2(:).device]);
YData = YData_([CIMR2(:).device]);
ZData = [CIMR2(:).value].';


%% further refine the edge color
XData = [XData.', 1.3, 1.7, -2.7, 2.5, 3.6, 3.7, 1.8].';
YData = [YData.', -0.8, -0.4, 0.8, 3.7, 2.6, 0.4, -0.9].';
ZData = [ZData.', 2, 2, 2, 2, 2, 2, 2].';
% Plot heat map
PlotHeatMap(XData,YData,ZData,1);
% Get the max
ZDataMax = max(ZData);

% Move graph to top
uistack(GraphFigure,'top');
%% Set Figure Lim
FigureMargin = 0.3;
% xmax = max(abs(XData));
% ymax = max(abs(YData));
% xymax = max(xmax,ymax);
% xlim([-xymax-FigureMargin,xymax+FigureMargin]);
% ylim([-xymax-FigureMargin,xymax+FigureMargin]);
xlim([min(GraphFigure.XData)-FigureMargin,max(GraphFigure.XData)+FigureMargin]);
ylim([min(GraphFigure.YData)-FigureMargin,max(GraphFigure.YData)+FigureMargin]);

%% Save
%print(figure(1),['Graph_' DataName '.png'],'-dpng','-r600');

colorbar('Ticks',[-2,-1,0,1,2],...
         'TickLabels',{'-2 very weak','-1 weak','0 normal','1 strong','2 very strong'})
