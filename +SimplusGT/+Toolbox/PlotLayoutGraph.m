% Plot the graph based on the nodal admittance matrix data
%
% Author(s): Yitong Li

function [Matrix,Data,Figure] = PlotLayoutGraph(Ybus)

% Convert Ybus to GraphMatrix
Matrix = SimplusGT.AbsMatrixElement(Ybus,'DiagFlag',0);

% Plot
Data = graph(Matrix,'upper');
Figure = plot(Data); grid on; hold on;

% Set
highlight(Figure,Data,'EdgeColor',[0,0,0],'LineWidth',1.1);       % Change all edges to black by default
highlight(Figure,Data,'NodeColor',[0,0,0]);                    	% Change all nodes to black by default
highlight(Figure,Data,'MarkerSize',4.5);
highlight(Figure,Data,'NodeFontSize',9);
highlight(Figure,Data,'NodeFontWeight','bold');

end