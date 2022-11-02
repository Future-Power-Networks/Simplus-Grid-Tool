% Author(s): Yitong Li

clear all
close all
clc

load_system('SimplusGT.slx');
Simulink.exportToVersion(bdroot,'Library\SimplusGT_2015a.slx','R2015a');
close_system('SimplusGT.slx');