clear all
clc
close all

MfilePath = mfilename('fullpath');
[RootPath,~,~]  = fileparts(MfilePath);
cd(RootPath);

addpath(genpath([RootPath,'/CustomerData']));	
addpath(genpath([RootPath,'/GenericFunction']));	
addpath(genpath([RootPath,'/Plot']));	
addpath(genpath([RootPath,'/SimulinkModel']));
addpath(genpath([RootPath,'/SystemObject']));
addpath(genpath([RootPath,'/Toolbox']));

fprintf('The toolbox is installed sucessfully! \n')
