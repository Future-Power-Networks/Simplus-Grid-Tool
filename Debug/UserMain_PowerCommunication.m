% This is the main function for analying the power grids by
% power-communication-isomorphism.

% Author(s): Yitong Li, Yunjie Gu

%% Prepare
clear all
clc
close all

%% Get color RGB
SimplusGT.ColorRgb();

%% Select data
UserData = 'NETS_NYPS_68Bus_Test';

%% Compare toolbox with nature
Enable_ComparisonToolbox = 0;           % 1/0: Compare the toolbox with nature
                                        % This enable value will also be used later when plotting
                                        % Notes:
                                        % The LPF of PLL will influence the comparison a lot
                                   
Enable_SmallSignalAnalysis = Enable_ComparisonToolbox;
SimplusGT.Toolbox.Main();

if Enable_ComparisonToolbox
    save('pole_sys.mat','pole_sys');
    clc
    close all
end

fprintf('\n')
fprintf('==================================\n')
fprintf('Synchronisation analysis.\n')
fprintf('==================================\n')

% Update power flow
[V,I] = SimplusGT.Communication.UpdateVI(PowerFlowNew);

%% Enable settings
% For testing participation analysis
Enable_FiedlerAbs               = 1;
Enable_NoneZeroXi               = 0;

% Enable control loop
Enable_VoltageNode_InnerLoop    = 1;    % 1/0: star-delta conversion for flux inductance of voltage node                ???
Enable_CurrentNode_InnerLoop    = 1;    % 1/0: inner-current loop impedance of current node                             ???

Enable_vq_PLL                   = 1;    % 1/0: change Q-PLL to vq-PLL
Enable_Change_Sign_PLL        	= 0;    % 1/0: change the sign of Q, epsilon_m = 90 or -90, for current node
Enable_PLL_LPF                  = 1;    % 1/0: if the PLL with an additional 100Hz LPF
w_tau                           = 2*pi*100;

Select_Ibus_Ref              	= 2;

% Enable plot
Enable_Plot_Eigenvalue          = 1;    % 1/0: Plot eigenvalues.
Enable_Plot_GraphOrigin         = 1;    % 1/0: Plot the graph of the original power system
Enable_Plot_GraphKH             = 0;    % 1/0: Plot the graph for KH
Enable_Plot_GraphAnalysis       = 0;    % 1/0: Plot the graph for analysis

% Initialize figure index
Fig_N = 2000;

%% Nodal admittance matrix
% Symbolic
s = sym('s');

% Calculate nodal admittance matrix
fprintf('Calculate nodal admittance matrix...\n')
W0 = Wbase;
Ybus = SimplusGT.Communication.YbusCalcSym(ListLineNew,W0,'albe');       % We get the alpha/beta frame admittance matrix

%% Reorder the Data

% Notes:
%
% In this code, the bus/node should be orderred in this sequence:
% [all voltage nodes, all current nodes, all floating bus nodes], i.e.,
% [v bus, ..., v bus, i bus, ..., i bus, f bus, ... f bus].
% Hence, in this section, we re-order the data obtained from excel first,
% to make sure that this required sequence can be obtained. Noting that,
% the device data, the power flow data, and the network line data should
% all be re-orderred.
%
% Maybe in the end of this code, I should also re-order the result back to
% its original sequence.

fprintf('Reorder the data...\n')
SimplusGT.Communication.ReorderData();

%% The Influence of Inner Loop on Ybus
fprintf('Evaluating the influence of bus type on node admittance matrix...\n')
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

% =============================
% Find the node index
% =============================

% Find the index of the first voltage bus/node
n_Vbus_1st = 1;

% Find the index of the first floating bus/node
n_Fbus_1st = N_Bus+1;  % Default: no floating bus
for i = 1:N_Bus
    if ApparatusSourceType(i) == 3
        n_Fbus_1st = i;
        break;
    end
end

% Check if the following buses are all floating buses
for i = n_Fbus_1st:N_Bus
    if ApparatusSourceType(i) ~= 3
        error(['Error']);
    end
end

% Find the index of the first current bus/node
n_Ibus_1st = N_Bus+1;   % Default: no current bus
for i = 1:N_Bus
    if ApparatusSourceType(i) == 2
        n_Ibus_1st = i;
        break;
    end
end

% Check if the following buses are all current buses
for i = n_Ibus_1st:(n_Fbus_1st-1)
    if ApparatusSourceType(i) ~= 2
        error(['Error']);
    end
end

% =============================
% Plot Graph
% =============================
if Enable_Plot_GraphOrigin
YbusOrigin = Ybus(Order_New2Old,Order_New2Old);
GraphMatrix = SimplusGT.Communication.NormMatrixElement(YbusOrigin,'DiagFlag',0);
Fig_N = Fig_N + 1;
figure(Fig_N)
GraphData = graph(GraphMatrix,'upper');
GraphFigure = plot(GraphData); grid on; hold on;
highlight(GraphFigure,GraphData,'EdgeColor','k','LineWidth',1);     % Change all edges nodes to black
highlight(GraphFigure,GraphData,'NodeColor','k');
highlight(GraphFigure,Index_Vbus,'NodeColor',RgbBlue);
highlight(GraphFigure,Index_Ibus,'NodeColor',RgbRed);
SaveGraphData{Fig_N} = GraphData;
SaveGraphFigure{Fig_N} = GraphFigure;
end

% =============================
% Update bus index
% =============================
if n_Ibus_1st>N_Bus
    if n_Fbus_1st <= N_Bus
        n_Ibus_1st = n_Fbus_1st;     % Update Ibus_1st
    end
end

% =============================
% Handle voltage node
% =============================
if Exist_Vbus == 0
    fprintf('Warning: The system has no voltage node.\n')
else
fprintf('Handling voltage node...\n')

for i = 1:(n_Ibus_1st-1)
J{i} = ApparatusParaNew{i}.J;
D{i} = ApparatusParaNew{i}.D;
J{i} = J{i}*Wbase;
D{i} = D{i}*Wbase;
% Notes: 
% Adding '*Wbase' is because the power swing equation rather than the
% torque swing equation is used, and P=T*w0 if w is not in per unit system.

Lsg = ApparatusParaNew{i}.wL;
Rsg = ApparatusParaNew{i}.R;
Zsg = s*Lsg + Rsg;
Y_sg{i} = 1/Zsg;

% Convert Ysg to double
Y_sg_{i} = double(subs(Y_sg{i},'s',1i*(W0+dW)));
Y_sg{i} = double(subs(Y_sg{i},'s',1i*W0));
end    

% Doing D-Y conversion
if Enable_VoltageNode_InnerLoop

% Prepare star-delta conversion by adding new buses
Ybus = SimplusGT.Communication.PrepareConvertDY(Ybus,n_Ibus_1st,N_Bus,Y_sg);
Ybus_ = SimplusGT.Communication.PrepareConvertDY(Ybus_,n_Ibus_1st,N_Bus,Y_sg_);
    
% Doing the star-delta conversion.
% Notes: Assume old voltage bus as zero current bus, and then switch the
% current and voltage for these buses so that current becomes input, and
% finally remove corresponding blocks because the input current is zero.
Ybus = SimplusGT.Communication.HybridMatrixYZ(Ybus,N_Bus+1);
Ybus_ = SimplusGT.Communication.HybridMatrixYZ(Ybus_,N_Bus+1);

% Eliminate the old voltage bus, i.e., zero current bus
Ybus = Ybus(1:N_Bus,1:N_Bus);
Ybus_ = Ybus_(1:N_Bus,1:N_Bus);

% Update V and I
% Notes: The steady-state voltage at voltage buses are changed if we split
% the inductor outside the apparatus. The "for loop" voltage calculation is
% equivalent to the matrix form "V = inv(Ybus)*I". But this matrix form
% would lead to wrong results in some cases when Ybus is almost
% non-invertible.
I = I(1:N_Bus,:);
for i = 1:(n_Ibus_1st-1)
    V(i) = V(i) + I(i)/Y_sg{i};
end

else
    fprintf('Warning: The voltage-node inner loop has been disabled.\n')
end

end

% =============================
% Handle current node
% =============================
if Exist_Ibus == 0
    fprintf('Warning: The system has no current node.\n')
else

fprintf('Handling current node...\n')

% Get the inner loop parameters

for i = n_Ibus_1st:(n_Fbus_1st-1)

% Notes: 
% All inverters have same current controllers
kp_i = ApparatusParaNew{i}.kp_i_dq;  
ki_i = ApparatusParaNew{i}.ki_i_dq;                                                    % ???
Lf   = ApparatusParaNew{i}.L;
Rf   = ApparatusParaNew{i}.R;

% Notes:
% The bandwidth of vq-PLL and Q-PLL would be different, because Q is
% proportional to id*vq. If we want to ensure the vq-PLL bandwidth of
% different inverters to be same, there Q-PLL bandwidth would probably be
% different because of different id output or active power output. We
% ensure the vq-PLL bandwidth here.
kp_pll{i} = ApparatusParaNew{i}.kp_pll;                                                  % ???
% kp_pll{i} = 0;
ki_pll{i} = ApparatusParaNew{i}.ki_pll;
PI_pll{i} = kp_pll{i} + ki_pll{i}/s;

end

wm = sym('wm');

% alpha/beta
Z_PIi = kp_i + ki_i/(s-1i*wm);        
Z_Lf = s*Lf+Rf;
Y_inv = (s-1i*wm)/((kp_i + s*Lf+Rf)*(s-1i*wm) + ki_i);

Y_inv = subs(Y_inv,'wm',W0);
Y_inv_ = double(subs(Y_inv,'s',1i*(W0+dW)));
Y_inv = double(subs(Y_inv,'s',1i*W0));

% Calculate Y_inv_prime
Y_inv_prime = (Y_inv_ - Y_inv)/(1i*dW);

% Notes:
% When current controller is very fast, Y_inv -> 0 and can be ignored.

% Add Y_inv to nodal admittance matrix
if Enable_CurrentNode_InnerLoop                                                 % ??? 
    
for i = n_Ibus_1st:(n_Fbus_1st - 1)
    % Self branch
    Ybus(i,i) = Ybus(i,i) + Y_inv;
    Ybus_(i,i) = Ybus_(i,i) + Y_inv_;
end

% Update I
I = Ybus*V;

else
    fprintf('Warning: The current-node inner loop has been disabled.\n')
end

end

% =============================
% Handle floating node
% =============================
% Notes:
% The floating bus (i.e., no device bus) is assumed as zero-current bus,
% and eliminated here after converting the Y matrix to Y-Z hybrid matrix.
if Exist_Fbus == 0
    fprintf('Warning: The system has no floating node.\n')
    YbusVIF = Ybus;
    YbusVI = Ybus;
else
    
fprintf('Eliminating floating node...\n')

YbusVIF = Ybus;

Ybus = SimplusGT.Communication.HybridMatrixYZ(Ybus,n_Fbus_1st);
Ybus_ = SimplusGT.Communication.HybridMatrixYZ(Ybus_,n_Fbus_1st);

N_Bus = n_Fbus_1st-1;
Ybus = Ybus(1:N_Bus,1:N_Bus);
Ybus_ = Ybus_(1:N_Bus,1:N_Bus);
YbusVI = Ybus;
V = V(1:N_Bus,:);
I = I(1:N_Bus,:);

end

%% Network matrix: K and Gamma
fprintf('Calculating network matrix: K and Gamma...\n')

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

% Update input and output to network admittance/impedance matrix, i.e.,
% Output = -Gbus*Input
Input = [V(1:n_Ibus_1st-1);
         I(n_Ibus_1st:end)];
Output = [I(1:n_Ibus_1st-1);
          V(n_Ibus_1st:end)];
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
ang_S_degree = angle(S)/pi*180;

% Get epsilon
for i = 1:N_Bus
    if ApparatusSourceType(i) == 1
        epsilon(i) = 0;         % W = P
    elseif ApparatusSourceType(i) == 2
        epsilon(i) = pi/2;      % W = Q
        
        if Enable_vq_PLL
            theta_i = angle(-I(i));
            theta_v = angle(V(i));
            epsilon(i) = pi/2 - (theta_i-theta_v);      % The Q direction is changed to vq direction.
        end
        
        if Enable_Change_Sign_PLL                                             
            if real( I(i)*exp(-1i*angle(V(i))) )>0                                                     
                epsilon(i) = -epsilon(i);
            end
            % Notes:
            % For current node, conventional PLL uses vq rather than Q to
            % achieve the synchronization. In load convention, Q = vq*id -
            % vd*iq. If iq = 0, Q = vq*id if iq = 0. This means, the power
            % flow id will influence the sign of Q, which also means the
            % sign of loop gain also has to be changed to ensure the system
            % stability. Here, we change the the value of epsilon_m
            % depending on the power flow, to ensure xi>=0 and the
            % stability. Noting that Q is in load convention, so, id should
            % also be load convention. That means, when I>0 which means the
            % IBR is generating active power, the sign of epsilon should be
            % changed.
        end
    else
        error(['Error']);
    end
end


% Get K matrix
for m = 1:N_Bus
    for n = 1:N_Bus
        if n ~= m
            ang_K(m,n) = epsilon(m) - angle(S(m,n)) - angle(Gbus(m,n));
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
[K11,K12,K21,K22] = SimplusGT.PartitionMatrix(K,n_Ibus_1st-1,n_Ibus_1st-1);
% Notes:
% For single-machine-load system, i.e., a single bus system, K = 0. How to
% understand this?

% Get Gamma matrix
for m = 1:N_Bus
    for n = 1:N_Bus
        ang_Gamma(m,n) = epsilon(m) - angle(S(m,n)) - angle(Gbus_prime(m,n));
        Gamma(m,n) = abs(S(m,n))*abs(Gbus_prime(m,n))*sin(ang_Gamma(m,n));
    end
end
Gamma = double(Gamma);
Gamma = -Gamma;    % For negative feedback
[Gamma11,Gamma12,Gamma21,Gamma22] = SimplusGT.PartitionMatrix(Gamma,n_Ibus_1st-1,n_Ibus_1st-1);

% Notes:
% K can be interpreted as the synchronizing torque coefficient, as dS is
% proportional to K*dtheta. Gamma can be interpreted as the damping torque
% coefficient, as dS is proportional to Gamma*dw.
%
% Pls be careful: The definition of K here is same to that in the paper.
% But the definition of gamma here is different from that in the paper.

%% Apparatus Matrix: T and H^{-1}
fprintf('Calculating apparatus matrix...\n')

% Initialize inertia matrix
Hinv = eye(N_Bus);      % Let Hinv be identity matrix initially

% Notes
% The representation of Hinv is very important, especially when considering
% the whole system KH. When seperately considerring KH_V and KH_I, H_V >>
% H_I, or equivalently, Hinv_V << Hinv_I, should be valid so that the
% voltage source is much slower than the current source.

% Choose the reference node
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
    
% Update inertia matrix for voltage node
for i = 1:(n_Ibus_1st-1)
    % The inertia of a SG is J
    Hinv(i,i) = 1/J{i};
    Hinv(i,i) = double(Hinv(i,i));
end
    
% Symbolic transfer function form:
% omega = 1/(D + J*s) * W;
% for i = 1:(n_Ibus_1st-1)
%     T_V_sym{i} = 1/(D{i}/J{i} + s);         % All T_V{i} should be same
% end
% F_V_sym = -s/T_V_sym{n_v_ref};

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
    
% Update Hinv for current node
for i = n_Ibus_1st:N_Bus
    % The inertia of an inverter is ki_pll.
    if Enable_PLL_LPF == 0                                                    	% ???
        Hinv(i,i) = ki_pll{i};
    else
        Hinv(i,i) = w_tau*kp_pll{i};
        ki_pll{i} = 0;
    end
    Hinv(i,i) = double(Hinv(i,i));
end
    
% Symbolic transfer function form
% omega = (kp_pll{i} + ki_pll{i}/s) * W
% for i = n_Ibus_1st:N_Bus
%     T_I_sym{i} = PI_pll{i}/ki_pll{i};       % All T_I{i} should be same
% end
% F_I_sym = -s/T_I_sym{n_i_ref};

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
tau = 1/w_tau;
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

%% K analysis
fprintf('K analysis...\n')

% KH
KH = Hinv*K;
[KH11,KH12,KH21,KH22] = SimplusGT.PartitionMatrix(KH,n_Ibus_1st-1,n_Ibus_1st-1);
[phi,xi,psi] = eig(KH);                 % phi is the right eigenvector matrix, psi is the left eigenvector matrix, xi is the eigenvalue matrix. Noting that phi^(-1) can also be regarded as a left eigenvector matrix.
phi_inv = phi^(-1);
xi_diag = diag(xi);
xi_diag_Hz = vpa(xi_diag/2/pi,5)
[xi_min,xi_min_index] = min( real(xi_diag) );
fprintf(['KH stability: ']);
if xi_min < -1e-5
    fprintf(['unstable xi.\n']);
    xi_min;
    xi_min_index;
else
    fprintf(['stable xi.\n']);
    % Converet x_min to 2nd smallest x_min, i.e., get the smallest non-zero
    % xi.
    if Enable_NoneZeroXi
        [xi_min2,xi_min2_index] = mink(real(xi_diag),2);
        xi_min = xi_min2(2);
        xi_min_index = xi_min2_index(2);
    end
end

% Power Flow Analysis
ListPowerFlowOrigin = ListPowerFlowNew(Order_New2Old,:);
PowerFlowAngle = ListPowerFlowOrigin(:,5);
[PowerFlowAngleMax,PowerFlowAngleMaxIndex] = max(abs(PowerFlowAngle));

% KH analysis
KH_r = KH(1,:);
KH_c = KH(:,1);
KH_Vec = KH_c.*(KH_r');
[KH_min,KH_min_Index] = min(abs(KH_Vec));

if Enable_Plot_GraphKH
GraphMatrix = NormMatrixElement(KH,'DiagFlag',0);
fig_n = fig_n + 1;
figure(fig_n)
PlotGraph();
SaveGraphData{fig_n} = GraphData;
SaveGraphFigure{fig_n} = GraphFigure;
end

% Print xi
xi_min_Hz = xi_min/2/pi
xi_min_index

% Right eigenvector
phi_r_xi = phi(:,xi_min_index);
phi_r_xi = phi_r_xi(Order_New2Old_NoFbus,1);
PhiRightMedian = median(phi_r_xi);
PhiRightMean = mean(phi_r_xi);
PhiRightPositive = find(phi_r_xi>=PhiRightMean);
PhiRightNegative = find(phi_r_xi<PhiRightMean);

if Enable_Plot_GraphAnalysis
GraphMatrix = NormMatrixElement(YbusOrigin,'DiagFlag',0);
Fig_N = Fig_N + 1;
figure(Fig_N)
PlotGraph();
SaveGraphData{Fig_N} = GraphData;
SaveGraphFigure{Fig_N} = GraphFigure;
highlight(GraphFigure,PhiRightPositive,'NodeColor',RgbYellow);
highlight(GraphFigure,PhiRightNegative,'NodeColor',RgbGreen);
end

% Left eigenvector
phi_l_xi = transpose(phi_inv(xi_min_index,:));
phi_l_xi = phi_l_xi(Order_New2Old_NoFbus,1);
PhiLeftPositive = find(phi_l_xi>=0);
PhiLeftNegative = find(phi_l_xi<0);

if Enable_Plot_GraphAnalysis
GraphMatrix = NormMatrixElement(YbusOrigin,'DiagFlag',0);
Fig_N = Fig_N + 1;
figure(Fig_N)
PlotGraph();
SaveGraphData{Fig_N} = GraphData;
SaveGraphFigure{Fig_N} = GraphFigure;
highlight(GraphFigure,PhiLeftPositive,'NodeColor',RgbYellow);
highlight(GraphFigure,PhiLeftNegative,'NodeColor',RgbGreen);
end

% Fiedler vector
% FiedlerVec = phi_r_xi.*phi_l_xi;
% FiedlerPositive = find(FiedlerVec>=0);
% FiedlerNegative = find(FiedlerVec<0);
% [~,FiedlerMinIndex] = min(FiedlerVec);
% [~,FiedlerMaxIndex] = max(abs(FiedlerVec))
FiedlerVec = phi_r_xi.*abs(phi_l_xi);
FiedlerMedian = median(FiedlerVec);
FiedlerAverage = mean(FiedlerVec);
FiedlerPositive = find(FiedlerVec>=FiedlerAverage);
FiedlerNegative = find(FiedlerVec<FiedlerAverage);
[~,FiedlerMinIndex] = min(FiedlerVec);
FiedlerAbsVec = abs(FiedlerVec);
[FiedlerAbsMax,FiedlerAbsMaxIndex] = max(FiedlerAbsVec);

if Enable_FiedlerAbs
    [FiedlerMaxN,FiedlerMaxNIndex] = maxk(FiedlerAbsVec,6);
else
    [FiedlerMaxN,FiedlerMaxNIndex] = maxk(FiedlerVec,6);
end

% Deal with the floating bus
for i = 1:length(FiedlerMaxNIndex)
    if FiedlerMaxNIndex(i) >= 50
        FiedlerMaxNIndex(i) = FiedlerMaxNIndex(i) + 17;
    elseif FiedlerMaxNIndex(i) >= 49
        FiedlerMaxNIndex(i) = FiedlerMaxNIndex(i) + 15;
    elseif FiedlerMaxNIndex(i) >= 46
        FiedlerMaxNIndex(i) = FiedlerMaxNIndex(i) + 13;
   	elseif FiedlerMaxNIndex(i) >= 44
        FiedlerMaxNIndex(i) = FiedlerMaxNIndex(i) + 11;
  	elseif FiedlerMaxNIndex(i) >= 34
      	FiedlerMaxNIndex(i) = FiedlerMaxNIndex(i) + 10;
   	elseif FiedlerMaxNIndex(i) >= 30
      	FiedlerMaxNIndex(i) = FiedlerMaxNIndex(i) + 9;
   	elseif FiedlerMaxNIndex(i) >= 29
      	FiedlerMaxNIndex(i) = FiedlerMaxNIndex(i) + 7;
   	elseif FiedlerMaxNIndex(i) >= 28
      	FiedlerMaxNIndex(i) = FiedlerMaxNIndex(i) + 5;
  	elseif FiedlerMaxNIndex(i) >= 21
      	FiedlerMaxNIndex(i) = FiedlerMaxNIndex(i) + 2;
   	elseif FiedlerMaxNIndex(i) >= 19
      	FiedlerMaxNIndex(i) = FiedlerMaxNIndex(i) + 1;
    end
end
FiedlerMaxN
FiedlerMaxNIndex

if Enable_Plot_GraphAnalysis
GraphMatrix = NormMatrixElement(YbusOrigin,'DiagFlag',0);
Fig_N = Fig_N + 1;
figure(Fig_N)
PlotGraph();
SaveGraphData{Fig_N} = GraphData;
SaveGraphFigure{Fig_N} = GraphFigure;
highlight(GraphFigure,FiedlerPositive,'NodeColor',RgbYellow);
highlight(GraphFigure,FiedlerNegative,'NodeColor',RgbGreen);
end

% GammaHphi
GammaHphi = inv(phi)*Hinv*Gamma*phi;
[GH11,GH12,GH21,GH22] = SimplusGT.PartitionMatrix(GammaHphi,n_Ibus_1st-1,n_Ibus_1st-1);
[~,sigma,~] = svd(GammaHphi);
sigma_max = max(max(sigma));

% Notes:
%
% We analyze the current node and voltage node seperately because their
% transfer functions T are different. By asuming the current node is much
% "faster" than the voltage node, we can analyze the current node first
% with assume no voltage node. Then, we can analyze the voltage node by
% adding current node back into K and gamma of the voltage node.
%
% If the power system is in the loop topology. KH and Gamma_Hphi should be
% semi-diagonalized and should have some relation. The values of KH will
% rotate in each row.
%
% If analyzing voltage and current nodes seperately, KH and Gamma_Hphi for
% two subloops should be re-calculated based on K, Hinv, Gamma, and can not
% be seperate directly from the whole system KH and Gamma_Hphi.

%% Calculate the state space representation
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
T12cl = feedback(T1cl,Gamma,feedin,feedout_L2);
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
fprintf('Plotting...\n')

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
legend('Yunjie Wihtout F Shift','Yunjie With F Shift','Toolbox')
end

%% Check stability
% By proposed criterion
% fprintf('Check the stability by the proposed criterion:\n')
% Stability Check Voltage Current
% StabilityCheckVoltageCurrent();

% By pole
fprintf('Check the stability by poles:\n')
UnstablePoleIndex = find(real(pole_T12cl)>1e-9);
UnstablePoleIndex0 = find(real(pole_T12cl)>0);
if isempty(UnstablePoleIndex)
    fprintf('Stable!\n')
else
    fprintf('Unstable!\n')
end
UnstablePole = pole_T12cl(UnstablePoleIndex)
RiskPole = pole_T12cl(UnstablePoleIndex0)

%% Save
% Admittance matrix
SaveData.YbusOrigin = YbusOrigin;
SaveData.YbusVI = YbusVI;
SaveData.YbusVIF = YbusVIF;

% Hybrid admittacne/impedance matrix
SaveData.GbusVI = GbusVI;
SaveData.GbusVIF = GbusVIF;

% Others
SaveData.KH = KH;
SaveData.FiedlerVec = FiedlerVec;
SaveData.Index_Vbus = Index_Vbus;
SaveData.Index_Ibus = Index_Ibus;
SaveData.Index_Fbus = Index_Fbus;
SaveData.Index_Ebus = Index_Ebus;
SaveData.Order_Old2New = Order_Old2New;
SaveData.Order_New2Old = Order_New2Old;
SaveData.pole_T1cl = pole_T1cl;
SaveData.pole_T12cl = pole_T12cl;

stop