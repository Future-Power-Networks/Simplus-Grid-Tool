% This is the main function for analying the power grids by
% communication theory

%% Prepare
clear all
clc
close all

%% Get color RGB
SimplusGT.ColorRgb();

%% Select data
% UserData = 'Test_68Bus_NETS_NYPS';      % Default NETS_NYPS system
% UserData = 'Test_68Bus_IBR_Load';       % IBRs with passvie loads
% UserData = 'Test_68Bus_IBR';            % IBRs with active loads
UserData = 'Test_68Bus_IBR_17';         % IBR at node 17 is repaced by a SG
% UserData = 'Test_68Bus_IBR_17_14';      % 17, 14
% UserData = 'Test_68Bus_IBR_17_14_7';    % 17, 14, 7

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
Enable_Plot_GraphOrigin         = 1;    % 1/0: Plot the graph of the original power system

% Initialize figure index
Fig_N = 2000;

%% Compare conventional small-signal (toolbox) with the proposed communication theory
Enable_ComparisonToolbox = 1;           % 1/0: Compare the toolbox with nature
                                        % This enable value will also be used later when plotting
                                        % Notes:
                                        % The LPF of PLL will influence the comparison a lot
             
%% Load data from excel by using toolbox functions
SimplusGT.Toolbox.Main();

if Enable_ComparisonToolbox
    save('pole_sys.mat','pole_sys');
end
clc
close all

%%
fprintf('\n')
fprintf('==================================\n')
fprintf('Synchronisation analysis.\n')
fprintf('==================================\n')

%% Update power flow
[V,I] = SimplusGT.Communication.UpdateVI(PowerFlowNew);

%% Nodal admittance matrix
% Symbolic
s = sym('s');

% Calculate nodal admittance matrix
fprintf('Calculate nodal admittance matrix...\n')
W0 = Wbase;
Ybus = SimplusGT.Communication.YbusCalcSym(ListLineNew,W0,'albe');       % We get the alpha/beta frame nodal admittance matrix

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

SimplusGT.Communication.ReorderData();

%% The Influence of node type and their parameters on Ybus
fprintf('Consider the influence of node type on node admittance matrix...\n')
% Nnotes:
% Ybus should statisfy: I = Ybus*V

% Convert Ybus to double
dW = 1e-10*(1+Wbase);
Ybus_ = double(subs(Ybus,'s',1i*(W0+dW)));   % Used for calculating derivative numerically
Ybus = double(subs(Ybus,'s',1i*W0));

% Notes:
% If using s-domain Ybus calculation and using vectors as input, then, the
% system admittance matrix has to be symmetric. Fortunately, the passive
% component, and the inner loops of inverters, are indeed symmtric in
% complex dq and alpha/beta frame.

% Plot Graph
if Enable_Plot_GraphOrigin
SimplusGT.Communication.PlotGraph();
end

% Find the node index
SimplusGT.Communication.FindNodeIndex();

% Handle voltage, current, and floating nodes
SimplusGT.Communication.HandleNode();

%% Network matrix
fprintf('Calculate network matrix: hybrid admittance/impedance matrix, or equivalently channel gain...\n')

% Convert the nodol admittance matrix to hybrid admittance/impedance matrix
Gbus = SimplusGT.Communication.HybridMatrixYZ(Ybus,n_Ibus_1st);
GbusVI  = Gbus;
GbusVIF = SimplusGT.Communication.HybridMatrixYZ(YbusVIF,n_Ibus_1st);

% For numerically calculating Gbus_prime later
Gbus_ = SimplusGT.Communication.HybridMatrixYZ(Ybus_,n_Ibus_1st);

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

%%
fprintf('Calculate network matrix: inertia, damping...\n')
% Initialize
Hmat = eye(N_Bus);
Hinv = inv(Hmat);
Dmat = eye(N_Bus);

% Update voltage node
if Exist_Vbus == 1
for i = 1:(n_Ibus_1st-1)
    % The inertia of a SG is J
    Hmat(i,i) = J{i};
    Dmat(i,i) = D{i};
    
  	Hinv(i,i) = 1/Hmat(i,i);
    Hinv(i,i) = double(Hinv(i,i));
end
end

% Update current node
if Exist_Ibus == 1
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

%%
fprintf('Calculating the linearized network matrix: K and FreqShift...\n')
% Get K matrix
for m = 1:N_Bus
    for n = 1:N_Bus
        if n ~= m
            ang_K(m,n) = mu(m) - angle(S(m,n)) - angle(Gbus(m,n));
            K(m,n) = abs(S(m,n))*abs(Gbus(m,n))*sin(ang_K(m,n));
        end
    end
end
for m = 1:N_Bus
    K_temp = 0;
    for n = 1:N_Bus
        if n~=m
            K_temp = K_temp - K(m,n);
        end
    end
    K(m,m) = K_temp;
end
K = double(K);
K = -K;                             % For negative feedback

% Get FreqShift matrix
for m = 1:N_Bus
    for n = 1:N_Bus
        ang_FreqShift(m,n) = mu(m) - angle(S(m,n)) - angle(Gbus_prime(m,n));
        FreqShift(m,n) = abs(S(m,n))*abs(Gbus_prime(m,n))*sin(ang_FreqShift(m,n));
    end
end
FreqShift = double(FreqShift);
FreqShift = -FreqShift;             % For negative feedback

% Notes:
% K can be interpreted as the synchronizing torque coefficient, as dS is
% proportional to K*dtheta. FreqShift can be interpreted as the damping
% torque coefficient, as dS is proportional to FreqShift*dw. K and
% FreqShift are all linearized results.

%% K analysis
fprintf('Analyze K...\n')
KH = Hinv*K;
[KH11,KH12,KH21,KH22] = SimplusGT.PartitionMatrix(KH,n_Ibus_1st-1,n_Ibus_1st-1);
[phi,xi,psi] = eig(KH);                 % phi is the right eigenvector matrix, psi is the left eigenvector matrix, xi is the eigenvalue matrix. Noting that phi^(-1) can also be regarded as a left eigenvector matrix.
phi_inv = phi^(-1);
xi_diag = diag(xi);
xi_diag = vpa(xi_diag,5)
[xi_min,xi_min_index] = min( real(xi_diag) );
fprintf(['KH stability: ']);
if xi_min < -1e-5
    fprintf(['unstable xi.\n']);
    xi_min;
    xi_min_index;
else
    fprintf(['stable xi.\n']);
    Enable_NoneZeroXi = 1;  % Converet x_min to 2nd smallest x_min, i.e., 
                            % get the smallest non-zero xi.
    if Enable_NoneZeroXi
        [xi_min2,xi_min2_index] = mink(real(xi_diag),2);
        xi_min = xi_min2(2);
        xi_min_index = xi_min2_index(2);
    end
end

if 0 
SimplusGT.Communication.AnalysisK();
end

%% Apparatus Matrix: T
fprintf('Calculate the apparatus state space model...\n')

% Notes
% The representation of Hinv is very important, especially when considering
% the whole system KH. When seperately considerring KH_V and KH_I, H_V >>
% H_I, or equivalently, Hinv_V << Hinv_I, should be valid so that the
% voltage source is much slower than the current source.

% Choose the reference node
Select_Ibus_Ref = 1;
if Select_Ibus_Ref == 1
    n_i_ref = n_Ibus_1st;  	% Select the first current node as the reference
elseif Select_Ibus_Ref == 2
    n_i_ref = N_Bus;        % Select the final current node as the reference
else
    error(['Error;']);
end
n_v_ref = 1;                % Select the first voltage node as the reference

% ================================
% Voltage node
% ================================
if Exist_Vbus == 1

% State space form:
% dx/dt = [domega]/dt = [-D/J,0]*[omega] + [1/J]*[W];
%         [dtheta]      [1   ,0] [theta]   [0  ]
% y = [omega] = [1,0]*[omega] + [0]*[W]
%     [theta]   [0,1] [theta]   [0]
Av = [-D{n_v_ref}/J{n_v_ref}, 0;
      1,                      0];
Bv = [1/J{n_v_ref};
      0];
Cv = [1,0;
      0,1];
Dv = [0;
      0];
T_V_ss = ss(Av,Bv,Cv,Dv);
T_V_ss = T_V_ss*J{n_v_ref};         % T_V_ss actually represents this system: 1/(D/J+s) without Hinv

end

% ================================
% Current node
% ================================
if Exist_Ibus == 1

% State space form
if Enable_PLL_LPF == 0                                                        	% ???
% Without PLL LPF
% dx/dt = [dw_pll_i]/dt = [0 ,0]*[w_pll_i] + [ki]*[W]
%         [dtheta  ]      [1, 0] [theta  ]   [kp]
% y = [omega] = [1, 0]*[w_pll_i] + [kp]*[W]
%     [theta]   [0 ,1] [theta  ]   [0 ]
Ai = [0, 0;
      1, 0];
Bi = [ki_pll{n_i_ref};
      kp_pll{n_i_ref}];
Ci = [1, 0;
      0, 1];
Di = [kp_pll{n_i_ref};
      0];
else
% With PLL LPF
% dx/dt = [dw      ]/dt
%         [dw_pll_i]
%         [dtheta  ]
% y = [omega]
%     [theta]
tau = 1/w_PLL_LPF;
Ai = [-1/tau,1/tau,0;
      0,0,0;
      1,0,0];
Bi = [kp_pll{n_i_ref}/tau;
      ki_pll{n_i_ref};
      0];
Ci = [1,0,0;
      0,0,1];
Di = [0;
      0];
end
T_I_ss = ss(Ai,Bi,Ci,Di);
T_I_ss = T_I_ss/Hinv(n_i_ref,n_i_ref);

end

%% Calculate the state space representation
fprintf('Calculate the whole system state space model...\n')
% Get whole system Tss
Tss = [[],[],[],[]];
for i = 1:N_Bus
    if ApparatusSourceType(i) == 1
        Tss = append(Tss,T_V_ss);
    elseif ApparatusSourceType(i) == 2
        Tss = append(Tss,T_I_ss);
    else
        error(['Error']);
    end
end

% Calculate the whole-system closed loop state space model
feedin = [1:N_Bus];
feedout_L1 = [1:N_Bus]*2;           % theta port
feedout_L2 = [1:N_Bus]*2-1;         % omega port
T1cl = feedback(Tss*Hinv,K,feedin,feedout_L1);
T12cl = feedback(T1cl,FreqShift,feedin,feedout_L2);
if ~isempty(T1cl.E) || ~isempty(T12cl.E)
    error(['Error: T1cl or T2cl is a dss system.']);
end
% T1cl = minreal(T1cl);
% T12cl = minreal(T12cl);

% Calculate the whole system pole
if 1
    [~,pole_T1cl] = eig(T1cl.A);
    pole_T1cl = diag(pole_T1cl)/2/pi;
    [~,pole_T12cl] = eig(T12cl.A);
    pole_T12cl = diag(pole_T12cl)/2/pi;
else
    pole_T1cl = pole(T1cl)/2/pi;
    pole_T12cl = pole(T12cl)/2/pi;
end

%% Plot
fprintf('Plot...\n')

% Plot: poles of state space system
if Enable_Plot_Eigenvalue
Fig_N = Fig_N+1;
figure(Fig_N)
scatter(real(pole_T1cl),imag(pole_T1cl),'x','LineWidth',1.5); hold on; grid on;
scatter(real(pole_T12cl),imag(pole_T12cl),'x','LineWidth',1.5); hold on; grid on;

legend('Loop1','Loop12')
end

if Enable_ComparisonToolbox
Fig_N = Fig_N+1;
figure(Fig_N)
scatter(real(pole_T1cl),imag(pole_T1cl),'x','LineWidth',1.5); hold on; grid on;
scatter(real(pole_T12cl),imag(pole_T12cl),'x','LineWidth',1.5); hold on; grid on;
pole_toolbox = load('pole_sys').pole_sys;
index = find(abs(imag(pole_toolbox))<35);
pole_toolbox = pole_toolbox(index);
index = find(real(pole_toolbox)>-1e3);
pole_toolbox = pole_toolbox(index);
scatter(real(pole_toolbox),imag(pole_toolbox),'x','LineWidth',1.5); hold on; grid on;
legend('Wihtout Freq Shift','With Freq Shift','Toolbox')
end

%% Check stability
fprintf('Check the stability by poles:\n')
UnstablePoleIndex     = find(real(pole_T12cl)>1e-9);
UnstablePoleIndexRisk = find(real(pole_T12cl)>0);
if isempty(UnstablePoleIndex)
    fprintf('Stable!\n')
else
    fprintf('Unstable!\n')
end
% UnstablePole = pole_T12cl(UnstablePoleIndex)
RiskPole = pole_T12cl(UnstablePoleIndexRisk)