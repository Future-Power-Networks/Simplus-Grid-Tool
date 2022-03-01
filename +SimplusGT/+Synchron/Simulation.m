% This function runs the simulation
%
% Author(s): Yitong Li

function [] = Simulation(Input,GAMMA,gamma,GAMMAFault,gammaFault,Hmat,Dmat,Wref,Wbase,N_Node,FigN)

options = odeset('MaxStep',1e-3);   % This para is very important and influences the results a lot.
x0_1 = angle(Input);
x0_2 = ones(N_Node,1)*Wbase;
x0 = [x0_1;
      x0_2];
  
StateSize = size(x0)
  
% Add a disturbance
% x0(1) = x0(1) - 5/180*pi;

switch 2
    case 1
        TimeRange = [0 10];
        [t_out,x_out] = ode45(@(t,x) SimplusGT.Synchron.StateEqu(x,GAMMA,gamma,Hmat,Dmat,Wref,Wbase),TimeRange,x0,options);
    case 2
        TimeRange = [0 10];
        [t_out,x_out] = ode45(@(t,x) SimplusGT.Synchron.StateEqu(x,GAMMA,gamma,Hmat,Dmat,Wref,Wbase),TimeRange,x0,options);
        t_out1 = t_out;
        x_out1 = x_out;
        clear('t_out','x_out');

        TimeRange = [10, 10.2];
        x0 = transpose(x_out1(end,:));
        [t_out,x_out] = ode45(@(t,x) SimplusGT.Synchron.StateEqu(x,GAMMAFault,gammaFault,Hmat,Dmat,Wref,Wbase),TimeRange,x0,options);
        t_out2 = t_out;
        x_out2 = x_out;
        clear('t_out','x_out');

        TimeRange = [10.2, 12];
        x0 = transpose(x_out2(end,:));
        [t_out,x_out] = ode45(@(t,x) SimplusGT.Synchron.StateEqu(x,GAMMA,gamma,Hmat,Dmat,Wref,Wbase),TimeRange,x0,options);
        t_out3 = t_out;
        x_out3 = x_out;
        clear('t_out','x_out');

        t_out = [t_out1;t_out2;t_out3];
        x_out = [x_out1;x_out2;x_out3];
    case 3
        TimeRange = [0, 8];
        [t_out,x_out] = ode45(@(t,x) SimplusGT.Synchron.StateEqu(x,GAMMA,gamma,Hmat,Dmat,Wref,Wbase),TimeRange,x0,options);
        t_out1 = t_out;
        x_out1 = x_out;
        clear('t_out','x_out');

        TimeRange = [8, 15];
        x0 = transpose(x_out1(end,:));
        x0(12) = x0(12) - pi/2;
        x0(28) = x0(28) - 0.05*Wbase;
        [t_out,x_out] = ode45(@(t,x) SimplusGT.Synchron.StateEqu(x,GAMMA,gamma,Hmat,Dmat,Wref,Wbase),TimeRange,x0,options);
        t_out2 = t_out;
        x_out2 = x_out;
        clear('t_out','x_out');

        t_out = [t_out1;t_out2];
        x_out = [x_out1;x_out2];
end

% Organize the data
theta_out = x_out(:,[1:N_Node]);
omega_out = x_out(:,[N_Node+1:2*N_Node]);
theta_out = theta_out - theta_out(:,1); % Get the angle difference
theta_out = theta_out/pi*180;
omega_out = omega_out/Wbase;

% theta_out = mod(theta_out,360);

TimeShift = 1.5;
TimeShift = 0;
TimeLimit = [0,2.5];

FigN = FigN + 1;
figure(FigN)
FigSize = [0.1 0.1 0.35 0.5];
set(gcf,'units','normalized','outerposition',FigSize);
subplot(2,1,1)
plot(t_out-TimeShift,theta_out);
% xlim(TimeLimit);
ylabel('Angle (Degree)')
subplot(2,1,2)
plot(t_out-TimeShift,omega_out);
% xlim(TimeLimit);
ylabel('Frequency (pu)')
xlabel('Time (s)')
% SimplusGT.mtit('Time-domain simulation');

% Notes:
% For time-domain simulation, Gbus should be positive and the convention
% should not be changed. Or the sign of Gamma should be changed.
%
% The time-domain simulation shows both the voltage and current
% angle/frequency rather than the voltage only.

end
