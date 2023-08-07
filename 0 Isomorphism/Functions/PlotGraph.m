% This function plots the graph for a power system

% Plot the graph
GraphData = graph(GraphMatrix,'upper');
GraphFigure = plot(GraphData); grid on; hold on;
highlight(GraphFigure,GraphData,'EdgeColor','k','LineWidth',1);     % Change all edges to black
highlight(GraphFigure,GraphData,'NodeColor','k');