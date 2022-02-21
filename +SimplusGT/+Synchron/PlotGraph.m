YbusOrigin = Ybus(Order_New2Old,Order_New2Old);
GraphMatrix = SimplusGT.Communication.NormMatrixElement(YbusOrigin,'DiagFlag',0);
Fig_N = Fig_N + 1;
figure(Fig_N)
GraphData = graph(GraphMatrix,'upper');
GraphFigure = plot(GraphData); grid on; hold on;
highlight(GraphFigure,GraphData,'EdgeColor','k','LineWidth',1); 	% Change all edges nodes to black
highlight(GraphFigure,GraphData,'NodeColor','k');
highlight(GraphFigure,Index_Vbus,'NodeColor',RgbBlue);              % Blue for voltage node
highlight(GraphFigure,Index_Ibus,'NodeColor',RgbRed);               % Red for current node
SaveGraphData{Fig_N} = GraphData;
SaveGraphFigure{Fig_N} = GraphFigure;