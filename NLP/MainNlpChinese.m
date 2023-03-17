clear all
clc
close all

% disp('Thanks for using Simplus Grid Tool! \n')
fprintf('=============================\n')
fprintf('欢迎使用Simplus Grid Tool！\n')
fprintf('=============================\n')

UserNlp.prompt = '请输入：';
WordLib.UserFunction = [{'小信号'},{'参与因子'},{'时域仿真'}];
WordLib.UserSystem = [{'默认'},{'IEEE'},{'14'},{'39'},{'57'},{''}]

fprintf('\n')
fprintf('请问您想使用的功能：\n')
fprintf('（全部功能；\n 小信号建模；参与因子；时域仿真...） \n')
fprintf('\n')
UserNlp.UserFunction = input(UserNlp.prompt,'s');
fprintf('\n')

fprintf('\n')
fprintf('请问您想测试的系统：\n')
fprintf('（默认系统UserData；单机无穷网系统；IEEE14节点；IEEE39节点；IEEE68节点...） \n')
fprintf('\n')
UserNlp.UserSystem = input(UserNlp.prompt,'s');
fprintf('\n')