% This is the main function for analying the power grids by
% communication theory

%% Prepare
clear all
clc
close all

%% Select data
% UserData = 'Test_68Bus_NETS_NYPS';      % Default NETS_NYPS system
% UserData = 'Test_68Bus_IBR_Load';       % IBRs with passvie loads
% UserData = 'Test_68Bus_IBR';            % IBRs with active loads
% UserData = 'Test_68Bus_IBR_17';         % IBR at node 17 is repaced by a SG
UserData = 'Test_68Bus_IBR_17_14';      % 17, 14
% UserData = 'Test_68Bus_IBR_17_14_7';    % 17, 14, 7

% UserData = 'Test_2Bus';
% UserData = 'Test_3Bus';

%% Enable settings
% Enable inner loop
Enable_VoltageNode_InnerLoop    = 1;    % 1/0: star-delta conversion for flux inductance of voltage node 
Enable_CurrentNode_InnerLoop    = 1;    % 1/0: inner-current loop impedance of current node

% PLL settings
Enable_vq_PLL                   = 1;    % 1/0: change Q-PLL to vq-PLL
Enable_PLL_LPF                  = 1;    % 1/0: if the PLL with an additional 100Hz LPF
w_PLL_LPF                       = 2*pi*100;     % Bandwidth of the PLL LPF

% Enable plot
Enable_Plot_Eigenvalue          = 1;    % 1/0: Plot eigenvalues.

% Initialize figure index
Fig_N = 2000;

%% Fault settings
NfaultBus = 3;
Tfault = 1/60*3;
             
%% Load data from excel by using toolbox functions
SimplusGT.Toolbox.Main();

%%
fprintf('\n')
fprintf('==================================\n')
fprintf('Synchronisation analysis.\n')
fprintf('==================================\n')

%% Update power flow
[V,I] = SimplusGT.Synchron.UpdateVI(PowerFlowNew);

%% Nodal admittance matrix
% Symbolic
s = sym('s');

% Calculate nodal admittance matrix
fprintf('Calculate nodal admittance matrix...\n')
Ybus = SimplusGT.Synchron.YbusCalcSym(ListLineNew,Wbase,'albe');       % We get the alpha/beta frame nodal admittance matrix

% Convert Ybus to double
dW = 1e-10*(1+Wbase);
Ybus_ = double(subs(Ybus,'s',1i*(Wbase+dW)));   % Used for calculating derivative numerically
Ybus = double(subs(Ybus,'s',1i*Wbase));
% Nnotes:
% Ybus should statisfy: I = Ybus*V

% Get Ybus during fault
ListLineFault = SimplusGT.Synchron.GetListLineFault(NfaultBus,ListLineNew,ListBus);
YbusFault = SimplusGT.Synchron.YbusCalcSym(ListLineFault,Wbase,'albe');
YbusFault = double(subs(YbusFault,'s',1i*Wbase));

%% Reorder the Data
fprintf('Reorder the data...\n')
% Notes:
%
% In this code, the bus/node should be orderred in this sequence:
% [all voltage nodes, all current nodes, all floating bus nodes], i.e.,
% [v bus, ..., v bus, i bus, ..., i bus, f bus, ... f bus].
% Hence, in this section, we re-order the data obtained from excel first,
% to make sure that this required sequence can be obtained. Noting that,
% the device data, the power flow data, and the network line data should
% all be re-orderred.

SimplusGT.Synchron.ReorderData();

%% The Influence of node type and their parameters on Ybus
fprintf('Consider the influence of node type on node admittance matrix...\n')

% Notes:
% If using s-domain Ybus calculation and using vectors as input, then, the
% system admittance matrix has to be symmetric. Fortunately, the passive
% component, and the inner loops of inverters, are indeed symmtric in
% complex dq and alpha/beta frame.

% Find the node index
SimplusGT.Synchron.FindNodeIndex();

% Handle voltage, current, and floating nodes
SimplusGT.Synchron.HandleNode();

%% Network matrix
fprintf('Calculate network matrix: hybrid admittance/impedance matrix, or equivalently channel gain...\n')

% Convert the nodol admittance matrix to hybrid admittance/impedance matrix
Gbus = SimplusGT.Synchron.HybridMatrixYZ(Ybus,n_Ibus_1st);
GbusVI  = Gbus;
GbusVIF = SimplusGT.Synchron.HybridMatrixYZ(YbusVIF,n_Ibus_1st);

% For numerically calculating Gbus_prime later
Gbus_ = SimplusGT.Synchron.HybridMatrixYZ(Ybus_,n_Ibus_1st);

% Notes:
% It should be ensured that the buses are listed in the form like this:
% [Vbus1, Vbus2, Vbus3, Ibus4, Ibus5, ...]
    
Gbus = -Gbus;  	% Change the power direction to load convention.
              	% Noting that this operation is different from Ybus
                % = -Ybus if the system has current nodes. The current
                % direction is not important actually.
ang_G_degree = angle(Gbus)/pi*180;
                
% For numerically calculating Gbus_prime
Gbus_ = -Gbus_;

% Get G_prime
% Notes: It is calculaed by numerical method
Gbus_prime = (Gbus_ - Gbus)/(1i*dW);         	% Consider

% Get fault Gbus
GbusFault = SimplusGT.Synchron.HybridMatrixYZ(YbusFault,n_Ibus_1st);
GbusFault = - GbusFault;

%% 
fprintf('Calculate network matrix: complex power...\n')
% Update input and output so that they correspond to the hybrid
% admittance/impedance matrix, i.e., Output = -Gbus*Input
Input = [V(1:n_Ibus_1st-1);
         I(n_Ibus_1st:end)];
Output = [I(1:n_Ibus_1st-1);
          V(n_Ibus_1st:end)];
      
% Normalize the current node because of PLL
InputNormalized = Input;        % Initialize
for i = 1:length(Input)
    if ApparatusSourceType(i) == 2
        % Notes:
        % If current source, then normalize the current of it, in order to
        % match the actual feedback signal, i.e., voltage rather than
        % power, of the PLL inverter later.
        %
        % This vq-PLL effect is considerred into the S matrix next, rather
        % than the T or H matrix. This effect only needs to be considerred
        % once.
        if Input(i) == 0
            InputNormalized(i,1) = 0;
        else
            InputNormalized(i,1) = Input(i)/abs(Input(i));
        end
    end
end

% Get S matrix
if Enable_vq_PLL
    S = conj(InputNormalized)*transpose(Input);
else
    S = conj(Input)*transpose(Input);
end

%% 
fprintf('Calculate network matrix: mu, GAMMA, and gamma...\n')
% Get mu
for i = 1:N_Bus
    if ApparatusSourceType(i) == 1          % Voltage node
        mu(i) = 0;         % W = P
    elseif ApparatusSourceType(i) == 2      % Current node
        mu(i) = pi/2;      % W = Q
        if Enable_vq_PLL
            theta_i = angle(-I(i));
            theta_v = angle(V(i));
            mu(i) = pi/2 - (theta_i-theta_v);      % The Q direction is changed to vq direction.
        end
    else
        error(['Error']);
    end
end

% Get GAMMA and gamma
for m = 1:N_Bus
    for n = 1:N_Bus
        GAMMA(m,n) = abs(Gbus(m,n)*S(m,n));
        gamma(m,n) = pi/2 + mu(m) - angle(Gbus(m,n));
    end
end

% Get GAMMA and gamma during fault
for m = 1:N_Bus
    for n = 1:N_Bus
        GAMMAFault(m,n) = abs(GbusFault(m,n)*S(m,n));
        gammaFault(m,n) = pi/2 + mu(m) - angle(GbusFault(m,n));
    end
end

%%
fprintf('Calculate network matrix: inertia, damping...\n')
% Initialize
Hmat = eye(N_Bus);
Hinv = inv(Hmat);
Dmat = eye(N_Bus);

% Update voltage node
if ExistVbus == 1
for i = 1:(n_Ibus_1st-1)
    % The inertia of a SG is J
    Hmat(i,i) = J{i};
    Dmat(i,i) = D{i};
    
  	Hinv(i,i) = 1/Hmat(i,i);
    Hinv(i,i) = double(Hinv(i,i));
end
end

% Update current node
if ExistIbus == 1
for i = n_Ibus_1st:N_Bus
    % The inertia of an inverter is ki_pll.
    if Enable_PLL_LPF == 0
        Hmat(i,i) = 1/ki_pll{i};                        % PI format
        Dmat(i,i) = kp_pll{i}/ki_pll{i};
    else                                    
        Hmat(i,i) = 1/(w_PLL_LPF*kp_pll{i});         	% LPF format
        Dmat(i,i) = 1/kp_pll{i};
        ki_pll{i} = 0;
    end
    Hinv(i,i) = 1/Hmat(i,i);
    Hinv(i,i) = double(Hinv(i,i));
end
end

% Power reference
for i = 1:N_Bus
    if ApparatusSourceType(i) == 1          % Voltage node
        Wref(i) = -PowerFlowNew{i}(1);
    elseif ApparatusSourceType(i) == 2      % Current node
        Wref(i) = 0;
 	end
end

%%
if 1
    SimplusGT.Synchron.SmallSignalAnalysis();
end

%%
if 0
fprintf('Run time-domain simulation...\n')
options = odeset('MaxStep',1e3);
TimeRange = [0 2];
x0_1 = angle(Input);
x0_2 = ones(N_Bus,1)*Wbase;
x0 = [x0_1;
      x0_2];
  
% Add a disturbance
% x0(1) = x0(1) - 5/180*pi;

if 1
    TimeRange = [0 4];
    [t_out,x_out] = ode45(@(t,x) SimplusGT.Synchron.StateEqu(x,GAMMA,gamma,Hmat,Dmat,Wref,Wbase),TimeRange,x0,options);
else
    [t_out,x_out] = ode45(@(t,x) SimplusGT.Synchron.StateEqu(x,GAMMA,gamma,Hmat,Dmat,Wref,Wbase),TimeRange,x0,options);
    t_out1 = t_out;
    x_out1 = x_out;
    clear('t_out','x_out');
    
    TimeRange = [2, 2.05];
    x0 = transpose(x_out1(end,:));
    [t_out,x_out] = ode45(@(t,x) SimplusGT.Synchron.StateEqu(x,GAMMAFault,gammaFault,Hmat,Dmat,Wref,Wbase),TimeRange,x0,options);
    t_out2 = t_out;
    x_out2 = x_out;
    clear('t_out','x_out');
    
  	TimeRange = [2.05, 4];
    x0 = transpose(x_out2(end,:));
    [t_out,x_out] = ode45(@(t,x) SimplusGT.Synchron.StateEqu(x,GAMMA,gamma,Hmat,Dmat,Wref,Wbase),TimeRange,x0,options);
    t_out3 = t_out;
    x_out3 = x_out;
    clear('t_out','x_out');
    
    t_out = [t_out1;t_out2;t_out3];
    x_out = [x_out1;x_out2;x_out3];
end

% Organize the data
theta_out = x_out(:,[1:N_Bus]);
omega_out = x_out(:,[N_Bus+1:2*N_Bus]);
theta_out = theta_out - theta_out(:,2); % Get the angle difference
theta_out = theta_out/pi*180;
omega_out = omega_out/Wbase;

% theta_out = mod(theta_out,360);

Fig_N = Fig_N + 1;
figure(Fig_N)
subplot(2,1,1)
plot(t_out,theta_out);
subplot(2,1,2)
plot(t_out,omega_out);

% Notes:
% For time-domain simulation, Gbus should be positive and the convention
% should not be changed. Or the sign of Gamma should be changed.
%
% The time-domain simulation shows both the voltage and current
% angle/frequency rather than the voltage only.

end