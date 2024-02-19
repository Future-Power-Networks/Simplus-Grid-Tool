clear all
clc
close all

load('NodeYZ.mat');
load('GraphFigure.mat');
x = GraphFigure.XData';
y = GraphFigure.YData';
v = NodeYZ;

fig = figure;

PlotHeatMap(x,y,v,fig,[0,0.6]);


