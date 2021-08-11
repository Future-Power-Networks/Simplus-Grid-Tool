# SimplusGT Manual for Users

## Introduction

This is the manual for toolbox users. If you have any inquires or face any problems when using the toolbox, please feel free to contact: Yitong Li (yitong.li15@imperial.ac.uk), Yunjie Gu (yg934@bath.ac.uk), and Yue Zhu (yue.zhu18@imperial.ac.uk).

## System Requirement

Matlab 2015a or later, with Simulink, Simscape/PowerSystem.

## Run Toolbox the First Time  

Run "InstallSimplusGT.m" by Matlab. That's all! The SimplusGT will automatically run and get results of an example power system saved in "UserData.xlsx". More examples can be found in "Examples" folder.

## Logic and Architecture of the Toolbox

The logic and architecture of this toolbox are shown in the figure below. Users only need to prepare a excel form "UserData.xlsx" and run "UserMain.m". Then, the tool will automatically do the static analysis, dynamic analysis, and time-domain simulation of the power system saved in the excel form.

![](https://raw.githubusercontent.com/Future-Power-Networks/Simplus-Grid-Tool/master/Documentations/Figures/Architecture.png)