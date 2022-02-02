% This is the main function for analying the power grids by
% communication theory

%% Prepare
clear all
clc
close all

%% Select data
% UserData = 'Test_68Bus_NETS_NYPS';        % Default NETS_NYPS system
% UserData = 'Test_68Bus_IBR';              % Full-IBR system
UserData = 'Test_68Bus_IBR_17_14';      	% IBRs at nodes 17, 14 are replaced by SGs
% UserData = 'Test_4Bus';                     % 4 bus system

%% Enable settings
% Enable inner loop
Enable_VoltageNode_InnerLoop    = 0;    % 1/0: star-delta conversion for flux inductance of voltage node 
Enable_CurrentNode_InnerLoop    = 0;    % 1/0: inner-current loop impedance of current node

% PLL settings
Enable_vq_PLL                   = 1;    % 1/0
Enable_PLL_LPF                  = 1;    % 1/0: if the PLL with an additional 100Hz LPF
w_PLL_LPF                       = 2*pi*100;     % Bandwidth of the PLL LPF

% Enable plot
Enable_Plot_Eigenvalue          = 1;    % 1/0: Plot eigenvalues.

% Initialize figure index
Fig_N = 2000;
             
%% Load data from excel by using toolbox functions
SimplusGT.Toolbox.Main();

%%
fprintf('\n')
fprintf('==================================\n')
fprintf('Synchronisation analysis.\n')
fprintf('==================================\n')

%% Nodal admittance matrix
% Symbolic
s = sym('s');

% Calculate nodal admittance matrix
fprintf('Calculate nodal admittance matrix...\n')
W0 = Wbase;
Ybus = SimplusGT.Communication.YbusCalcSym(ListLineNew,W0,'albe');       % We get the alpha/beta frame nodal admittance matrix

% Convert Ybus to numerical value
dW = 1e-10*(1+Wbase);
Ybus_ = double(subs(Ybus,'s',1i*(W0+dW)));   % Used for calculating derivative numerically
Ybus = double(subs(Ybus,'s',1i*W0));

% Get V and I
[V,I] = SimplusGT.Communication.UpdateVI(PowerFlowNew);

% Nnotes:
% Ybus should statisfy: I = Ybus*V

%% Reorder the Data
fprintf('Reorder the data...\n')
% Notes:
%
% In this code, the bus/node should be orderred in this sequence:
% [all voltage nodes, all current nodes, all floating/empty bus nodes], i.e.,
% [v bus, ..., v bus, i bus, ..., i bus, f bus, ... f bus].
% Hence, in this section, we re-order the data obtained from excel first,
% to make sure that this required sequence can be obtained. Noting that,
% the device data, the power flow data, and the network line data should
% all be re-orderred.

SimplusGT.Communication.ReorderData();

%% The Influence of node type and their parameters on Ybus
fprintf('Consider the influence of node type on node admittance matrix...\n')

% Notes:
% If using s-domain Ybus calculation and using vectors as input, then, the
% system admittance matrix has to be symmetric. Fortunately, the passive
% component, and the inner loops of inverters, are indeed symmtric in
% complex dq and alpha/beta frame.

% Find the node index
SimplusGT.Communication.FindNodeIndex();

% Handle voltage, current, and floating/empty nodes
SimplusGT.Communication.HandleNode();

%% Network matrix
fprintf('Calculate network matrix: hybrid admittance/impedance matrix, or equivalently channel gain...\n')

% Convert the nodol admittance matrix to hybrid admittance/impedance matrix
Gbus = SimplusGT.Communication.HybridMatrixYZ(Ybus,n_Ibus_1st);
Gbus_ = SimplusGT.Communication.HybridMatrixYZ(Ybus_,n_Ibus_1st);

Gbus = -Gbus;  	% Change the power direction to load convention.
              	% Noting that this operation is different from Ybus
                % = -Ybus if the system has current nodes. The current
                % direction is not important actually.

Gbus_prime = (Gbus_ - Gbus)/(1i*dW);                % Consider 

%% 
fprintf('Calculate network matrix: complex power...\n')
% Update input and output so that they correspond to the hybrid
% admittance/impedance matrix, i.e., Output = Gbus*Input
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
      
%% 
fprintf('Calculate network matrix: mu, GAMMA, and gamma...\n')
% Get mu
for i = 1:N_Bus
    if ApparatusSourceType(i) == 1          % Voltage node
        mu(i,:) = 0;         % W = P
    elseif ApparatusSourceType(i) == 2      % Current node
        % mu(i,:) = -pi/2;      % W = -Q
        mu(i,:) = pi/2;
        % mu(i,:) = 0;
        
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

% Output
% mu
% GAMMA
% gamma

%%
fprintf('Calculate network matrix: inertia, damping...\n')
% Initialize
Hmat = eye(N_Bus);
Dmat = eye(N_Bus);
Hinv = inv(Hmat);

% Update voltage node
if Exist_Vbus == 1
for i = 1:(n_Ibus_1st-1)
    % The inertia of a SG is J
    Hmat(i,i) = J{i};
    Dmat(i,i) = D{i};
    
    Hinv(i,i) = 1/Hmat(i,i);
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
end
end

%% 
SimplusGT.Communication.SmallSignalAnalysis();